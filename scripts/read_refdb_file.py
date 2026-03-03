import csv
import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple
import numpy as np

FITS = ['bayes','theilsen','ransac','quantile','tukey']
RCS = ['avg','luk','pou','sch','wan','wis']
LACS_PATH = '/Users/kumaranbaskaran/lacs/'
def read_refdb_file(filename):
    data = []
    with open(filename) as f:
        reader = csv.reader(f,delimiter='\t')
        for row in reader:
            data.append(row)
    print (len(data))
    out= 'queue bmrbid in'
    print (data[0])
    refdb_ha=data[0].index('HA_CSDIFF')
    refdb_ca=data[0].index('CA_CSDIFF')
    refdb_cb=data[0].index('CB_CSDIFF')
    refdb_c=data[0].index('CO_CSDIFF')
    refdb_n= data[0].index('N_CSDIFF')
    for row in data[1:]:
        #print (len(row),row[0],row[0][3:])
        #out=f'{out} {row[0][3:]}'
        bmrb_id = row[0][3:]
        refdb_ca_offset = row[refdb_ca]
        refdb_cb_offset = row[refdb_cb]
        refdb_c_offset = row[refdb_c]
        refdb_n_offset = row[refdb_n]
        refdb_ha_offset = row[refdb_ha]
        for fit in FITS:
            for rc in RCS:
                lacs_output = Path(f'{LACS_PATH}{fit}/{rc}/{bmrb_id}/{bmrb_id}_{fit}.json')
                print (read_json_dict(lacs_output,fit))

def parse_json_file(path: Path) -> Tuple[Path, Dict[str, Any] | None, str | None]:
    """
    Try to load a JSON file and return (path, data, error_message).
    On success, error_message is None. On failure, data is None.
    """
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return path, data, None
    except json.JSONDecodeError as e:
        return path, None, f"JSON decode error at line {e.lineno} col {e.colno}: {e.msg}"
    except OSError as e:
        return path, None, f"OS error: {e.strerror} ({e.filename})"
    except Exception as e:
        return path, None, f"Unexpected error: {type(e).__name__}: {e}"



def read_json_dict(path: Path, fit) -> Dict[str, Any]:
    p,data,err = parse_json_file(path)
    print (err,p)
    if data is not None:
        print (data)
        for list_id in data:
            if fit in ['theilsen','ransac','quantile','tukey']:
                for atom in data[list_id][f'offsets']:
                    print(atom, data[list_id][f'offsets'][atom])
            else:
                for atom in data[list_id][f'offsets_{fit}']:
                    print (atom,data[list_id][f'offsets_{fit}'][atom])

if __name__ == '__main__':
    read_refdb_file('bmrpdbnew.txt')