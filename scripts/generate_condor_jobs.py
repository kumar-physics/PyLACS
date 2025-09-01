



def generate_condor_submit_file(fname):
    f=open(fname,'w')
    f.write("executable \t= /home/nmrbox/kbaskaran/SecondaryStructureStatistics.py\n")
    f.write("arguments \t= $(filename)\n")
    f.write("error \t= err.$(Process)\n")
    f.write("output \t= out.$(Process)\n")
    f.write("log = log.$(Process)\n")
    pair_list = _get_bmrb_pdb_mapping()
    flist=""
    for k in pair_list:
        bmrb = k['bmrb_id']
        for pdb in k['pdb_ids']:
            flist+=f'{pdb.lower()}-{bmrb} '
    f.write(f'queue filename in {flist}\n')
    f.close()