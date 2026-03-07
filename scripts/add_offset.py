import sys

import pynmrstar
import os
import json
import numpy as np
FITS = ['bayes','theilsen','ransac','quantile','tukey']
RCS = ['avg','luk','pou','sch','wan','wis']
def add_offset(nmrstar_file,atom,offset):
    ent = pynmrstar.Entry.from_file(nmrstar_file)
    for saveframe in ent:
        if saveframe.category == 'assigned_chemical_shifts':
            for loop in saveframe.loops:
                if loop.category == '_Atom_chem_shift':
                    colums = loop.tags
                    #print (colums)
                    atom_idx = colums.index('Atom_ID')
                    val_idx = colums.index('Val')
                    for row in range(len(loop.data)):
                        if loop.data[row][atom_idx] in ['CA','CB','C','N']:
                            loop.data[row][val_idx]= float(loop.data[row][val_idx])-offset
    fo=open('/Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test/test.str','w')
    fo.write(str(ent))
    return True

def run_rnage(nmrstar_file):
    fo=open('controlled_benchmark2.csv','w')
    for offset in np.arange(3.5,4.0,0.5):
        for fit in FITS:
            for rc in RCS:
                x=add_offset(nmrstar_file,'CA',offset)
                if rc == 'avg':
                    cmd = f'source /Users/kumaranbaskaran/Projects/bmrb/PyLACS/venv/bin/activate \n /Users/kumaranbaskaran/Projects/bmrb/PyLACS/pylacs/src/pylacs/lacs.py /Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test/test.str --data-id test --out /Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test --method={fit}'
                else:
                    cmd = f'source /Users/kumaranbaskaran/Projects/bmrb/PyLACS/venv/bin/activate \n /Users/kumaranbaskaran/Projects/bmrb/PyLACS/pylacs/src/pylacs/lacs.py /Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test/test.str --data-id test --out /Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test --method={fit} --rc-model={rc}'

                #os.system('source /Users/kumaranbaskaran/Projects/bmrb/PyLACS/venv/bin/activate')
                os.system(cmd)
                with open(f'/Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/test/test_{fit}.json','r') as f:
                    if fit != 'bayes':
                        data = json.load(f)
                        ca_off = data['1']['offsets']['ca']
                        cb_off = data['1']['offsets']['cb']
                        c_off = data['1']['offsets']['c']
                        n_off = data['1']['offsets']['n']
                    else:
                        data = json.load(f)
                        ca_off = data['1']['offsets_bayes']['ca']['mean']
                        cb_off = data['1']['offsets_bayes']['cb']['mean']
                        c_off = data['1']['offsets_bayes']['c']['mean']
                        n_off = data['1']['offsets_bayes']['n']['mean']

                fo.write(f'{fit},{rc},{offset},{ca_off},{cb_off},{c_off},{n_off}\n')
    fo.close()

if __name__ == '__main__':
    nmrstar_file = '/Users/kumaranbaskaran/Projects/bmrb/PyLACS/scratch/bmr19710_3.str'
    #add_offset(nmrstar_file,'CA',2.3)
    run_rnage(nmrstar_file)