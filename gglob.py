from optparse import OptionParser
import json
import glob

from sys import stderr
import os
from wnd_cp_data import wnd_cp_indiv
import time
import numpy as np

class gglob:
    """
    a class for loading and saving giant matrices of DTS 
    windowed genome data. A gglob_dir has 2 components:
    1. gglob.index
    2. chr{*}.npz for all chrs 
    """
        
    def __init__(self, gglob_dir, contig):
        #open up the index
        idx_data = None
        
        with open("%s/gglob.idx"%gglob_dir) as F:
            idx_data = json.load(F)
        
        self.indivs = idx_data['indivs'] 
        self.wnd_size = idx_data['wnd_size'] 
        self.wnd_slide = idx_data['wnd_slide'] 
        self.contig = contig
       
        self.wnd_starts = np.load("%s/%s.wnd_starts.npy"%(o.gglob_dir, contig))
        self.wnd_wnds = np.load("%s/%s.wnd_ends.npy"%(o.gglob_dir, contig))
        self.cp_matrix = np.load("%s/%s.cp_matrix.npy"%(o.gglob_dir, contig))
        
        self.sunk_wnd_starts = np.load("%s/%s.sunk_wnd_starts.npy"%(o.gglob_dir, contig))
        self.sunk_wnd_wnds = np.load("%s/%s.sunk_wnd_ends.npy"%(o.gglob_dir, contig))
        self.sunk_cp_matrix = np.load("%s/%s.sunk_cp_matrix.npy"%(o.gglob_dir, contig))
        
        #print "%s/%s.npz"%(gglob_dir, contig)
        #npz = np.load("%s/%s.npz"%(gglob_dir, contig))
        #print npz.files
        #print npz['wnd_starts'] 
        #print npz['wnd_ends'] 
        #print npz['cp_matrix'] 

        #self.wnd_starts = npz['wnd_starts']
        #self.wnd_ends = npz['wnd_ends']
        #self.cp_matrix = npz['cp_matrix']
        #
        #self.sunk_wnd_starts = npz['sunk_wnd_starts']
        #self.sunk_wnd_ends = npz['sunk_wnd_ends']
        #self.sunk_cp_matrix = npz['sunk_cp_matrix']

        assert self.sunk_cp_matrix.shape[0] == len(self.indivs) 
            
if __name__=="__main__":
        
    opts = OptionParser()

    usage_txt = "python gglob.py [options] [--init || --setup_chr {chr}]"
    opts = OptionParser(usage=usage_txt)

    opts.add_option('','--in_DTS_dir',dest='DTS_dir', default=None)
    opts.add_option('','--in_sunk_DTS_dir',dest='sunk_DTS_dir', default=None)

    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--sunk_contigs',dest='fn_sunk_contigs', default=None)
    
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    
    opts.add_option('','--gglob_dir',dest='gglob_dir')

    opts.add_option('','--setup_chr',dest='setup_chr', default=None)
    opts.add_option('','--init',dest='init', action='store_true', default=False)


    (o, args) = opts.parse_args()

    if o.init == False and o.setup_chr == None:
        opts.print_help()

    elif o.init:
        indivs = []
        idx_data = {} 
        for f in glob.glob("%s/%s*"%(o.DTS_dir,o.DTS_prefix)):
            indiv = f.split("/")[-1].replace(o.DTS_prefix,"") 
            if os.path.exists("%s/%s%s"%(o.sunk_DTS_dir, o.DTS_prefix, indiv )):
                indivs.append(indiv)
            else:
                print >>stderr, "skiping: %s - no associated sunk DTS"
        
        idx_data = {"indivs":indivs, "wnd_size":o.wnd_size, "wnd_slide":o.wnd_slide}

        with open("%s/gglob.idx"%o.gglob_dir, 'w') as F:
            json.dump(idx_data, F)
    
    else:
        
        idx_data = None
        with open("%s/gglob.idx"%o.gglob_dir) as F:
            idx_data = json.load(F)
        
        DTS_pre="%s/%s"%(o.DTS_dir, o.DTS_prefix) 
        sunk_DTS_pre="%s/%s"%(o.sunk_DTS_dir, o.DTS_prefix) 
        
        wnd_size = idx_data['wnd_size']
        indivs = idx_data['indivs']
        n_indivs = len(indivs)
        contig = o.setup_chr
        
        t = time.time()
        rand_wnd_cp = wnd_cp_indiv("%s%s"%(DTS_pre, indivs[0]), o.fn_contigs, wnd_size)
        wnd_starts, wnd_ends = rand_wnd_cp.get_wnds_by_chr(contig)
        cp_matrix = np.zeros((n_indivs, wnd_starts.shape[0]))

        rand_sunk_wnd_cp = wnd_cp_indiv("%s%s"%(sunk_DTS_pre, indivs[0]), o.fn_sunk_contigs, wnd_size)
        sunk_wnd_starts, sunk_wnd_ends = rand_sunk_wnd_cp.get_wnds_by_chr(contig)
        sunk_cp_matrix = np.zeros((n_indivs, sunk_wnd_starts.shape[0]))
        
        correct = not (contig in ["chrY", "chrX"])

        for i, indiv in enumerate(indivs):
            
            wnd_cp = wnd_cp_indiv("%s%s"%(DTS_pre, indiv),
                                  o.fn_contigs,
                                  wnd_size)
            
            cp_matrix[i,:] = wnd_cp.get_cps_by_chr(contig, correct=correct) 

            sunk_wnd_cp = wnd_cp_indiv("%s%s"%(sunk_DTS_pre, indiv), 
                                      o.fn_sunk_contigs,
                                      wnd_size)
            
            sunk_cp_matrix[i,:] = sunk_wnd_cp.get_cps_by_chr(contig, correct=correct) 
        
        np.save("%s/%s.wnd_starts"%(o.gglob_dir, contig), wnd_starts)
        np.save("%s/%s.wnd_ends"%(o.gglob_dir, contig), wnd_ends)
        np.save("%s/%s.cp_matrix"%(o.gglob_dir, contig), cp_matrix)
        
        np.save("%s/%s.sunk_wnd_starts"%(o.gglob_dir, contig), sunk_wnd_starts)
        np.save("%s/%s.sunk_wnd_ends"%(o.gglob_dir, contig), sunk_wnd_ends)
        np.save("%s/%s.sunk_cp_matrix"%(o.gglob_dir, contig), sunk_cp_matrix)
        

        #np.savez_compressed("%s/%s"%(o.gglob_dir, contig), wnd_starts=wnd_starts,
        #                                                   wnd_ends=wnd_ends, 
        #                                                   cp_matrix=cp_matrix,
        #                                                   sunk_wnd_starts=sunk_wnd_starts,
        #                                                   sunk_wnd_ends=sunk_wnd_ends,
        #                                                   sunk_cp_matrix=sunk_cp_matrix)
        print >>stderr, "done (%f)"%(time.time()-t)


