from optparse import OptionParser
import glob
import pandas as pd
import gzip


"""
as far as I can tell this program uses 33.5g of RAM? 
that's a-lot. tooo bloody much maybe I can just zcat?  
"""

if __name__=="__main__":
    

    print "use merge_individual_tables.sh instead, it is faster and more memory efficient"
    opts = OptionParser()
    opts.add_option('', '--in_dir', dest='in_dir')
    opts.add_option('', '--outfile', dest='fn_out')
    (o, args) = opts .parse_args()
    
    all_calls =  []
    header = ['rank', 'chr', 'start', 'end', 'mu', 'p', 'adjusted_p', 'window_size']   
    tlen = len(glob.glob("%s/*.gz"%(o.in_dir)))
    for k, f in enumerate(glob.glob("%s/*.gz"%(o.in_dir))):
        calls = pd.read_csv(f,sep = '\t', header=False, compression = 'gzip', names = header )
        all_calls.append(calls)
        fname = f.split("/dCGH_")[-1]
        g_test, g_ref = fname.split(".")[0:2]
        l = len(calls['rank'])
        calls['indiv_reference'] =  pd.Series([g_ref for i in xrange(l)], index = calls.index)
        calls['indiv_test'] =  pd.Series([g_test for i in xrange(l)], index = calls.index)
        print "%d/%d (%f%%) %s"%(k, tlen, 100*float(k)/tlen,  f)
     
    call_table = pd.concat(all_calls) 
    #call_table.to_csv(open(o.fn_out),
    #F_gz = gzip.open(o.fn_out,'w')
    #F_out = open(o.fn_out,'w')
    #call_table.to_csv(F_out,sep='\t',index=False)
    #F_gz.close()
    #F_out.close()
    call_table.to_hdf(o.fn_out, 'table', append=False, complib='zlib', mode='w', complevel=9, fletch32=True)

    

        
