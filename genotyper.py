import glob


from wnd_cp_data import wnd_cp_indiv


class genotyper:

    def __init__(self, 
                 dts_dir, 
                 fn_contigs, 
                 wnd_size, 
                 F_fnToIndiv=lambda x: x.split("/")[-1].replace("500_bp_","")):
        
        self.wnd_cp_by_indiv={}

        for fn_DTS in glob.glob("%s/*"%o.dts_dir):

            indiv = F_fnToIndiv(fn_DTS)
            self.wnd_cp_by_indiv[indiv] = wnd_cp_indiv(fn_DTS,
                                                       fn_contigs,
                                                       wnd_size)




