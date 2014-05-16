

class info_io(object):

    def __init__(self, FOUT):
        self.FOUT = FOUT
        self.fields = ['contig','start','end']
        self.fields+= ["mu_mu_d", 
                       "max_mu_d", 
                       "min_mu_d", 
                       "f_correct_direction",
                       "min_z",
                       "wnd_size", 
                       "bic_delta",
                       "n_clusts",
                       "min_allele_count",
                       "Lscore",
                       "min_inter_label_dist",
                       "singleton_P"]
    
        outstr =  "\t".join(self.fields)
        self.FOUT.write("%s\n"%outstr)

    def init_entry(self):
        new_entry = {}
        return {}

    def update_entry(self, entry, key,value):
        entry[key] = value 
    
    def update_entries(self, entry, key_vals):
        for k, v in key_vals.iteritems():
            entry[k] = v
    
    def output_entry(self, entry):
        outstr = "\t".join(["{0}".format(entry[f]) for f in self.fields])
        self.FOUT.write("%s\n"%outstr)
