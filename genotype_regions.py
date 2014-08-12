#!/usr/bin/env python
"""
A simple utility for getting copy number from each DenseTrackSet or gglob for each region in a given file.
"""
import argparse
from collections import OrderedDict
import csv
import os
from wnd_cp_data import wnd_cp_indiv
import json
import pandas as pd
import numpy as np
import pdb

def get_regions(regions_file, supported_contigs):
    '''Load regions by chromosome with starts and ends for each region in
    separate lists and a placeholder for copy number results by sample name.'''

    with open(regions_file, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
                                              
        regions_by_chromosome = OrderedDict()
        for row in reader:                  
            if row[0] not in supported_contigs:
                continue                      
                                             
            if row[0] not in regions_by_chromosome: 
                regions_by_chromosome[row[0]] = {  
                    "starts": [],                 
                    "ends": []                   
                }                               
                                               
            regions_by_chromosome[row[0]]["starts"].append(int(row[1]))
            regions_by_chromosome[row[0]]["ends"].append(int(row[2]))
    return regions_by_chromosome

def gglob_get_cp_by_regions(cps, wnd_starts, wnd_ends, contig, start, end):
    """
    Get median copy number for each region. Based on genotype function
    from gglob_genotyper
    """

    l = cps.shape[1]
    wnd_start = np.searchsorted(wnd_starts[0], start)
    wnd_end = min(np.searchsorted(wnd_ends[0], end)+1,l-1)
    if wnd_end<=wnd_start:
        wnd_end=wnd_start+1
    copies = np.median(cps.ix[:, wnd_start:wnd_end+1],1)
    return copies.tolist()

def genotype_gglob_regions(gglob_dir, regions_file, contigs, window, sunk):
    """
    Get copy number for each individual in each region in the gglob dir for WSSD or SUNK.
    """
    if contigs is None:
        supported_contigs = ['chr' + str(x) for x in range(1,23) + ['X', 'Y']]
    else:
        with open(contigs, 'r') as fh:
            supported_contigs = set([line.strip().split('\t')[0] for line in fh])

    regions_by_chromosome = get_regions(regions_file, supported_contigs)

    idx_data = None

    with open("%s/gglob.idx" % gglob_dir) as F:
        idx_data = json.load(F)
        
    sample_names = idx_data['indivs']
    wnd_size = idx_data['wnd_size']
    wnd_slide = idx_data['wnd_slide']

    yield ["chromosome", "start", "end"] + sample_names

    sunk_string = 'sunk_' if sunk else ''

    for chromosome, regions in regions_by_chromosome.iteritems():
        cp_matrix = pd.read_hdf('%s/%s.%s%s.h5' % (gglob_dir, chromosome, sunk_string, 'cp_matrix'), '%scp_matrix' % sunk_string)
        wnd_ends = pd.read_hdf('%s/%s.%s%s.h5' % (gglob_dir, chromosome, sunk_string, 'wnd_ends'), '%swnd_ends' % sunk_string)
        wnd_starts = pd.read_hdf('%s/%s.%s%s.h5' % (gglob_dir, chromosome, sunk_string, 'wnd_starts'), '%swnd_starts' % sunk_string)

        print "Genotyping %s" % chromosome
        for i in xrange(len(regions['starts'])):
            copies = gglob_get_cp_by_regions(cp_matrix, wnd_starts, wnd_ends, chromosome, 
                                             regions['starts'][i], regions['ends'][i])
            if copies is not None:
                yield list([chromosome, str(regions['starts'][i]), str(regions['ends'][i])] + 
                       copies)
            else:
                yield []

def genotype_DTS_regions(dts_list_file, regions_file, contigs, window):
    """
    Get copy number from each DenseTrackSet for each region in the given file.
    """
    # Load supported contigs.
    with open(contigs, "r") as fh:
        supported_contigs = set([line.strip().split("\t")[0] for line in fh])

    regions_by_chromosome = get_regions(regions_file, supported_contigs)


    # Load DTS file list.
    with open(dts_list_file, "r") as fh:
        dts_list = [line.strip() for line in fh]

    # Build a sorted list of sample names.
    dts_list.sort()
    sample_names = [os.path.basename(dts) for dts in dts_list]

    yield ["chromosome", "start", "end"] + sample_names

    for chromosome, regions in regions_by_chromosome.iteritems():
        copies_by_sample = {}

        for i in xrange(len(dts_list)):
            dts = dts_list[i]
            sample_name = sample_names[i]
            sample = wnd_cp_indiv(dts, contigs, window)

            copies_by_sample[sample_name] = sample.get_cp_by_regions(
                chromosome,
                regions["starts"],
                regions["ends"]
            )

        for i in xrange(len(regions["starts"])):
            yield list([chromosome, str(regions["starts"][i]), str(regions["ends"][i])] +
                       [str(copies_by_sample[sample][i]) for sample in sample_names])


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--gglob_dir", help="directory of gglob files")
    parser.add_argument("--dts_list_file", help="list of paths to DTS files to genotype")
    parser.add_argument("regions_file", help="BED file of regions to genotype")
    parser.add_argument("genotypes_file", help="output file for genotypes")
    parser.add_argument("--contigs", default="/net/eichler/vol4/home/jlhudd/projects/DTS_caller/500_bp_windows.pkl.contigs")
    parser.add_argument("--window", type=int, default=500)
    parser.add_argument("--sunk", action='store_true')
    args = parser.parse_args()

    genotypes_prefix = args.genotypes_file.split('.')[0]

    if args.gglob_dir is None and args.dts_list_file is None:
        print "Must specify gglob_dir or dts_list_file.\n"
        exit(1)

    if args.gglob_dir is not None:
        genotype_func = genotype_gglob_regions
    else:
        genotype_func = genotype_DTS_regions

    with open(args.genotypes_file, "w") as oh:
        writer = csv.writer(oh, delimiter="\t")
        if args.dts_list_file is not None:
            for row in genotype_DTS_regions(args.dts_list_file, args.regions_file, args.contigs, args.window):
                writer.writerow(row)
        elif args.gglob_dir is not None: 
            for row in genotype_gglob_regions(args.gglob_dir, args.regions_file, args.contigs, args.window, sunk=args.sunk):
                writer.writerow(row)
                
