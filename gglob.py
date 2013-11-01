from optparse import OptionParser
import json
import glob

from sys import stderr
import os
from wnd_cp_data import wnd_cp_indiv
import time
import numpy as np
import pandas as pd

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
       
        keys = ["wnd_starts","wnd_ends","cp_matrix","sunk_wnd_starts","sunk_wnd_ends","sunk_cp_matrix"]
        fn_in = "%s/%s"%(gglob_dir,contig)

        mats_by_key = {} 
        for k in keys:
            stderr.write("loading %s..."%k)
            stderr.flush()
            t=time.time()
            df = pd.read_hdf("%s.%s.h5"%(fn_in,k),k)
            mats_by_key[k] = df.as_matrix()
            stderr.write("done (%fs)\n"%(time.time()-t))

        self.wnd_starts = mats_by_key['wnd_starts'][:,0]
        self.wnd_ends = mats_by_key['wnd_ends'][:,0]
        self.cp_matrix = mats_by_key['cp_matrix']
        
        self.sunk_wnd_starts = mats_by_key['sunk_wnd_starts'][:,0]
        self.sunk_wnd_ends = mats_by_key['sunk_wnd_ends'][:,0]
        self.sunk_cp_matrix = mats_by_key['sunk_cp_matrix']
        
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
        
        fn_out =  "%s/%s"%(o.gglob_dir, contig)
        
        mats = {"wnd_starts":wnd_starts,
                "wnd_ends":wnd_ends,
                "cp_matrix":cp_matrix,
                "sunk_wnd_starts":sunk_wnd_starts,
                "sunk_wnd_ends":sunk_wnd_ends,
                "sunk_cp_matrix":sunk_cp_matrix
                }
         
        for k, mat in mats.iteritems():
            print >>stderr, "writing out %s..."%k
            t=time.time()
            df = pd.DataFrame(mat)
            df.to_hdf("%s.%s.h5"%(fn_out,k),k,complevel=1,complib='zlib')
            #df.to_hdf("%s.%s.h5"%(fn_out,k),k,complevel=5,complib='lzo')
            #df.to_hdf("%s.%s.h5"%(fn_out,k),k)
            print >>stderr, "done (%fs)"%(time.time()-t)
        
        #np.save("%s/%s.wnd_starts"%(o.gglob_dir, contig), wnd_starts)
        #np.save("%s/%s.wnd_ends"%(o.gglob_dir, contig), wnd_ends)
        #np.save("%s/%s.cp_matrix"%(o.gglob_dir, contig), cp_matrix)
        #
        #np.save("%s/%s.sunk_wnd_starts"%(o.gglob_dir, contig), sunk_wnd_starts)
        #np.save("%s/%s.sunk_wnd_ends"%(o.gglob_dir, contig), sunk_wnd_ends)
        #np.save("%s/%s.sunk_cp_matrix"%(o.gglob_dir, contig), sunk_cp_matrix)
        

        #np.savez_compressed("%s/%s"%(o.gglob_dir, contig), wnd_starts=wnd_starts,
        #                                                   wnd_ends=wnd_ends, 
        #                                                   cp_matrix=cp_matrix,
        #                                                   sunk_wnd_starts=sunk_wnd_starts,
        #                                                   sunk_wnd_ends=sunk_wnd_ends,
        #                                                   sunk_cp_matrix=sunk_cp_matrix)
        print >>stderr, "done (%f)"%(time.time()-t)


