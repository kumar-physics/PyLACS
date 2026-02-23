import re
import re
import os
from pathlib import Path
import csv
import plotly.express as px
import pandas as pd

def parse_reference_offsets(file_path):
    """
    Extract all atom offsets (CO, CA, CB, N, etc.) from 'Detected reference offsets' section.
    Returns a dict {atom_name: offset_value}.
    """
    offsets = {}
    capture = False
    buffer = []

    header_pattern = re.compile(r"Detected\s+reference\s+offsets", re.IGNORECASE)
    atom_pattern = re.compile(r"\b([A-Z][A-Z0-9]*)\s*[:=]?\s*([-+]?\d*\.\d+|\d+)", re.IGNORECASE)

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # Start capturing after header
            if header_pattern.search(line):
                capture = True
                continue

            # Stop when section ends
            if capture:
                if not line.strip():
                    break
                if re.match(r"^\s*[A-Z].*:", line) and not atom_pattern.search(line):
                    break
                buffer.append(line.strip())

    if not buffer:
        return {}

    text = " ".join(buffer)
    matches = atom_pattern.findall(text)
    for atom, val in matches:
        try:
            offsets[atom] = float(val)
        except ValueError:
            continue

    return offsets


def parse_all_out_files(folder_path, recursive=False):
    """
    Parse all .out files in the given folder and return {filename: {atom: offset}}
    """
    results = {}
    folder = Path(folder_path)

    pattern = "**/*.out" if recursive else "*.out"
    out_files = list(folder.glob(pattern))

    if not out_files:
        print(f"No .out files found in {folder_path}")
        return results

    for file in out_files:
        offsets = parse_reference_offsets(file)
        if offsets:
            results[file.name] = offsets

    # Print summary
    print("\n=== Parsed Reference Offsets Summary ===\n")
    panav_results={}
    for fname, atoms in results.items():
        #print(f"{fname.split('.')[0]}:")
        for atom, val in atoms.items():
            #print(f" {fname.split('.')[0]},{atom},{val}")
            if fname.split('.')[0] not in panav_results:
                panav_results[fname.split('.')[0]]={}
            if atom == 'CA':
                panav_results[fname.split('.')[0]]['ca']=val
            elif atom == 'CB':
                panav_results[fname.split('.')[0]]['cb'] = val
            elif atom == 'CO':
                panav_results[fname.split('.')[0]]['c'] = val
            elif atom == 'N':
                panav_results[fname.split('.')[0]]['n'] = val
            else:
                print (fname.split('.')[0],atom,val)
        #print()
    print (len(panav_results))
    lacs_results = read_lacs_results('lacs_results_1026.csv')
    print (len(lacs_results))
    missing=[]
    panav={'ca':[],'cb':[],'c':[],'n':[]}
    lacs={'ca':[],'cb':[],'c':[],'n':[]}
    id = []

    for bid in panav_results:
        if bid not in lacs_results:
            missing.append(bid)
        else:
            id.append(bid)
            for atom in ['ca','cb','c','n']:
                try:
                    pval=float(panav_results[bid][atom])
                except KeyError:
                    pval = None
                except ValueError:
                    pval = None
                try:
                    lval=float(lacs_results[bid][atom])
                except KeyError:
                    lval = None
                except ValueError:
                    lval = None
                panav[atom].append(pval)
                lacs[atom].append(lval)

    # fo=open('PAVAV_LACS_offsets.csv','w')
    # for i in range(len(id)):
    #     fo.write(f'{id[i]},{lacs['ca'][i]},{lacs['cb'][i]},{lacs['c'][i]},{lacs['n'][i]},{panav['ca'][i]},{panav['cb'][i]},{panav['c'][i]},{panav['n'][i]}\n')
    # fo.close()
    flg = ['Same' if abs(lacs['ca'][i]-panav['ca'][i])<=0.5 else 'Diff' for i in range(len(lacs['ca']))]
    s=flg.count('Same')
    d=flg.count('Diff')
    df = pd.DataFrame({
        'BMRB_ID':id,
        'LACS (CA)' : lacs['ca'],
        'LACS (CB)': lacs['cb'],
        'LACS (C)': lacs['c'],
        'LACS (N)': lacs['n'],
        'PANAV (CA)': panav['ca'],
        'PANAV (CB)': panav['cb'],
        'PANAV (C)': panav['c'],
        'PANAV (N)': panav['n'],
        'Diff': [f'Same{s}' if abs(lacs['ca'][i]-panav['ca'][i])<=0.5 else f'Diff{d}' for i in range(len(lacs['ca']))]
    })
    fig = px.scatter(df,x='LACS (CA)',y= 'PANAV (CA)',hover_name='BMRB_ID',color='Diff')
    fig.write_html('PANAV_vs_LACS_CA.html')
    fig = px.scatter(df, x='LACS (CB)', y='PANAV (CB)', hover_name='BMRB_ID')
    fig.write_html('PANAV_vs_LACS_CB.html')
    fig = px.scatter(df, x='LACS (C)', y='PANAV (C)', hover_name='BMRB_ID')
    fig.write_html('PANAV_vs_LACS_C.html')
    fig = px.scatter(df, x='LACS (N)', y='PANAV (N)', hover_name='BMRB_ID')
    fig.write_html('PANAV_vs_LACS_N.html')
    fig = px.scatter(df, x='PANAV (CA)', y='PANAV (CB)', hover_name='BMRB_ID')
    fig.write_html('PANAV_CA_CB.html')
    fig = px.scatter(df, x='LACS (CA)', y='LACS (CB)', hover_name='BMRB_ID')
    fig.write_html('LACS_CA_CB.html')
    # fig = px.histogram(df,x='Diff')
    # fig.show()
    print (len(sorted(missing)))
    return panav_results


def read_lacs_results(fname):
    with open(fname,'r') as f:
        csv_reader = csv.reader(f)
        lacs_result = {}
        for row in csv_reader:
            if row[0] not in lacs_result:
                lacs_result[row[0]]={}
            lacs_result[row[0]]={'ca':row[1],'cb':row[2],'c':row[3],'n':row[4]}
    return lacs_result

#
# def parse_reference_offsets(file_path):
#     """
#     Extract all atom offsets (CO, CA, CB, N, etc.) from 'Detected reference offsets' section.
#     Works with multi-line tables and various spacings.
#     """
#     offsets = {}
#     capture = False
#     buffer = []
#
#     header_pattern = re.compile(r"Detected\s+reference\s+offsets", re.IGNORECASE)
#     atom_pattern = re.compile(r"\b([A-Z][A-Z0-9]*)\s*[:=]?\s*([-+]?\d*\.\d+|\d+)", re.IGNORECASE)
#
#     with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
#         for line in f:
#             # Start capturing after header
#             if header_pattern.search(line):
#                 capture = True
#                 continue
#
#             # Stop when section ends
#             if capture:
#                 if not line.strip():
#                     # blank line → section ended
#                     break
#                 if re.match(r"^\s*[A-Z].*:", line) and not atom_pattern.search(line):
#                     # a new header like "Summary:" or "File:" etc.
#                     break
#                 buffer.append(line.strip())
#
#     if not buffer:
#         print("No 'Detected reference offsets' section found or it’s empty.")
#         return {}
#
#     # Join all lines and extract all atom:value pairs
#     text = " ".join(buffer)
#     matches = atom_pattern.findall(text)
#     for atom, val in matches:
#         try:
#             offsets[atom] = float(val)
#         except ValueError:
#             continue
#
#     # Debug: show what was captured
#     print("Captured lines under 'Detected reference offsets':")
#     print("\n".join(buffer))
#     print("\nExtracted Offsets:\n")
#
#     for atom, val in offsets.items():
#         print(f"{atom}: {val}")
#
#     return offsets
#


if __name__ == "__main__":
    folder = "/Users/kumaranbaskaran/panav_output/split_output"  # <-- change path as needed
    parse_all_out_files(folder, recursive=False)
    # read_lacs_results('lacs_results_1026.csv')


# if __name__ == "__main__":
#     file_path = '/Users/kumaranbaskaran/panav_output/split_output/4998_1.out'
#     results = parse_reference_offsets(file_path)
#     for atom, offset in results.items():
#         print(f"{atom}: {offset}")
