import os
import json

import fontTools.t1Lib
import pandas as pd
import plotly.express as px
import pynmrstar
from typing import Dict, List, Optional, Sequence, Tuple

def read_str_files_in_a_folder(root_dir):
    """
        Reads STR files inside every folder under the given root directory.

        Parameters:
            root_dir (str): Path to the root directory
        """
    before = {}
    after = {}
    id = []
    atms = ['ca', 'cb', 'c', 'n']
    fo_bf = open('/home/nmrbox/kbaskaran/entropy/cs_before.csv','w')
    fo_af = open('/home/nmrbox/kbaskaran/entropy/cs_after.csv', 'w')
    for atom in atms:
        before[atom] = []
        after[atom] = []
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            bmrb_id = file.split("_")[0]
            rebox_path = f'/reboxitory/2025/08/BMRB/macromolecules/bmr{bmrb_id}/bmr{bmrb_id}_3.str'
            if file.lower().endswith(".str"):
                file_path = os.path.join(dirpath, file)
                print (file_path,file, bmrb_id)
                bf=read_star(rebox_path)
                af=read_star(file_path)
                for k in bf:
                    fo_bf.write(f"{','.join(k)},{bf[k]}\n")
                for k in af:
                    fo_af.write(f"{','.join(k)},{af[k]}\n")

ResidueKey = Tuple[str, str, str, str,str,str,str]
def read_star(file_name: str) -> Dict[ResidueKey, float]:
    """Read an NMR-STAR file and extract atom chemical shifts.

    Parameters
    ----------
    file_name : str
        Path to a ``.str`` file.

    Returns
    -------
    dict
        Mapping ``{list_id -> {(entity, assembly, comp_index, comp_id) -> {atom -> value}}}``.

    Notes
    -----
    Only rows in the ``Atom_chem_shift`` category are parsed. Non-numeric
    values are skipped.
    """
    try:
        ent = pynmrstar.Entry.from_file(file_name)
    except FileNotFoundError:
        return {}

    loops = ent.get_loops_by_category('Atom_chem_shift')
    if not loops:
        return {}
    col = loops[0].get_tag_names()

    idx_entity = col.index('_Atom_chem_shift.Entity_ID')
    idx_assy = col.index('_Atom_chem_shift.Entity_assembly_ID')
    idx_comp_index = col.index('_Atom_chem_shift.Comp_index_ID')
    idx_comp = col.index('_Atom_chem_shift.Comp_ID')
    idx_atom = col.index('_Atom_chem_shift.Atom_ID')
    idx_val = col.index('_Atom_chem_shift.Val')
    idx_list = col.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')

    cs_data: Dict[ResidueKey, float] = {}
    for cs_loop in loops:
        if not cs_loop.data:
            continue
        list_id = cs_loop.data[0][idx_list]
        #cs_data.setdefault(list_id, {})
        for row in cs_loop.data:
            key: ResidueKey = (ent.entry_id,row[idx_entity], row[idx_assy], row[idx_comp_index], row[idx_comp],list_id,row[idx_atom])
            try:
                if row[idx_atom] in ['CA','CB','C','N']:
                    cs_data[key] = float(row[idx_val])
            except (TypeError, ValueError):
                continue
    return cs_data

def read_json_in_folders(root_dir):
    """
    Reads JSON files inside every folder under the given root directory.

    Parameters:
        root_dir (str): Path to the root directory
    """
    before={}
    after={}
    id=[]
    atms=['ca','cb','c','n']
    for atom in atms:
        before[atom]=[]
        after[atom]=[]
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.lower().endswith(".json"):
                file_path = os.path.join(dirpath, file)
                try:
                    with open(file_path, "r") as f:
                        (head, tail) = os.path.split(file_path)
                        corrected_file = f'{os.path.split(head)[0]}/corrected/{tail.split("_")[0]}/{tail}'
                        if os.path.exists(corrected_file):
                            data = json.load(f)
                            print(f"\n--- Contents of {os.path.basename(file_path)} ---")
                            (head,tail)=os.path.split(file_path)
                            print (head,tail,os.path.split(head))
                            #json_data = json.dumps(data)
                            #print(json.dumps(data, indent=4))
                            for list_id in data:
                                id.append(tail.split("_")[0])
                                for atom in atms:
                                    try:
                                        before[atom].append(data[list_id]['offsets'][atom])
                                    except KeyError:
                                        before[atom].append(None)
                            corrected_file = f'{os.path.split(head)[0]}/corrected/{tail.split("_")[0]}/{tail}'
                            print (corrected_file)
                            with open(corrected_file,'r') as f2:
                                data2 = json.load(f2)
                                for list_id in data2:
                                    for atom in atms:
                                        try:
                                            after[atom].append(data2[list_id]['offsets'][atom])
                                        except KeyError:
                                            after[atom].append(None)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

    print (len(id))
    for atom in atms:
        print (len(before[atom]),len(after[atom]))
    df = pd.DataFrame({
        "BMRB_ID" : id,
        "CA(before)":before['ca'],
        "CA(after)":after['ca'],
        "CB(before)": before['cb'],
        "CB(after)": after['cb'],
        "C(before)": before['c'],
        "C(after)": after['c'],
        "N(before)": before['n'],
        "N(after)": after['n']
    })
    fig = px.scatter(df,x='CA(before)',y="CA(after)",hover_name="BMRB_ID",marginal_x="histogram",marginal_y="histogram")
    fig.show()
    fig.write_html('ca_correlation.html')
    fig = px.scatter(df, x='CB(before)', y="CB(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('cb_correlation.html')
    fig = px.scatter(df, x='C(before)', y="C(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('c_correlation.html')
    fig = px.scatter(df, x='N(before)', y="N(after)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('n_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="CB(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_cb_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="C(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_c_correlation.html')
    fig = px.scatter(df, x='CA(before)', y="N(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('ca_n_correlation.html')
    fig = px.scatter(df, x='C(before)', y="N(before)", hover_name="BMRB_ID", marginal_x="histogram",
                     marginal_y="histogram")
    fig.show()
    fig.write_html('c_n_correlation.html')

# Example usage:
if __name__ == "__main__":
    root_directory = "/home/nmrbox/kbaskaran/lacs/tukey"
    read_json_in_folders(root_directory)
