

class info_io(object):

    def __init__(self, FOUT):
        self.FOUT = FOUT
        self.fields = ['contig','start','end']
        self.header_init = False
    
    def init_header(self, entry):
        
        for k,v in entry.iteritems():
            if not k in self.fields:
                self.fields.append(k)

        outstr =  "\t".join(self.fields)
        self.FOUT.write("%s\n"%outstr)

        self.header_init = True

    def init_entry(self):
        new_entry = {}
        return {}

    def update_entry(self, entry, key,value):
        entry[key] = value 
    
    def update_entries(self, entry, key_vals):
        for k, v in key_vals.iteritems():
            entry[k] = v
    
    def output_entry(self, entry):

        if not self.header_init:
            self.init_header(entry)
        
        outstr = "\t".join(["{0}".format(entry[f]) for f in self.fields])
        self.FOUT.write("%s\n"%outstr)



