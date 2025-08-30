import csv
import json
from pylacs.lacs import run_lacs
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def read_csv_file(csvfile):
    flg=False
    with open(csvfile,'r') as csvFile:
        csvF = csv.reader(csvFile)
        for row in csvF:
            bmrb_id = row[0]
            print (bmrb_id)
            if bmrb_id == '5168':
                flg = True
            if flg and bmrb_id not in ['4077','24152','5168']:
                str_file = f'/reboxitory/2025/08/BMRB/macromolecules/bmr{bmrb_id}/bmr{bmrb_id}_3.str'
                data_id = bmrb_id
                for method in ["tukey","theilsen","ransac","quantile","bayes"]:
                    out_dir = Path(f'../scratch/{data_id}')
                    plots = True
                    report=run_lacs(str_file=str_file,method=method,data_id=data_id,outdir=out_dir,plots=plots)
                    with open(out_dir /f'{data_id}_{method}.json', "w", encoding="utf-8") as f:
                        json.dump(report, f, indent=2)


if __name__ == "__main__":
    read_csv_file('Table_S1_clean.csv')