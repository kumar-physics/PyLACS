#!/usr/bin/env python3
"""
This script is python version og LACS software  https://link.springer.com/article/10.1007/s10858-005-1717-0
"""
import pynmrstar
import plotly.express as px
import plotly.graph_objects as go
from random_coil import RandomCoil
import numpy as np

ONE_TO_THREE = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                               'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                               'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                               'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
THREE_TO_ONE = dict([(value, key) for key, value in ONE_TO_THREE.items()])


def read_str_file(file_name):
    ent = pynmrstar.Entry.from_file(file_name)
    chemical_shift_loops = ent.get_loops_by_category('Atom_chem_shift')
    column_names = chemical_shift_loops[0].get_tag_names()
    cs_data ={}
    entity_idx = column_names.index('_Atom_chem_shift.Entity_ID')
    entity_assembly_idx = column_names.index('_Atom_chem_shift.Entity_assembly_ID')
    comp_index_idx = column_names.index('_Atom_chem_shift.Comp_index_ID')
    comp_idx = column_names.index('_Atom_chem_shift.Comp_ID')
    atom_idx = column_names.index('_Atom_chem_shift.Atom_ID')
    val_idx = column_names.index('_Atom_chem_shift.Val')
    list_idx = column_names.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')
    for cs_loop in chemical_shift_loops:
        if cs_loop.data[0][list_idx] not in cs_data:
            cs_data[cs_loop.data[0][list_idx]]={}
        for row in cs_loop.data:
            if (row[entity_idx],row[entity_assembly_idx],
                              row[comp_index_idx],row[comp_idx]) not in cs_data[cs_loop.data[0][list_idx]]:
                cs_data[cs_loop.data[0][list_idx]][(row[entity_idx],row[entity_assembly_idx],
                              row[comp_index_idx],row[comp_idx])] ={}
            cs_data[cs_loop.data[0][list_idx]][(row[entity_idx], row[entity_assembly_idx],
                                                row[comp_index_idx], row[comp_idx])][row[atom_idx]]=float(row[val_idx])

    return cs_data

def lacs(str_file,rc_model = None):
    cs_data = read_str_file(str_file)
    rc_shifts = RandomCoil()
    atom_list = rc_shifts.atoms()
    lacs_data = {}
    for cs_list in cs_data:
        if cs_list not in lacs_data:
            lacs_data[cs_list] = {'d_ca':[],'d_cb':[],'d_ca_cb':[],'tag':[]}
        for residue in cs_data[cs_list]:
            if residue[-1] in THREE_TO_ONE:
                try:
                    delta_ca = cs_data[cs_list][residue]['CA']-rc_shifts.get_value(residue[-1],'CA',rc_name=rc_model)
                    delta_cb = cs_data[cs_list][residue]['CB']-rc_shifts.get_value(residue[-1],'CB',rc_name=rc_model)
                    # delta_c = cs_data[cs_list][residue]['C']-rc_shifts.get_value(residue[-1],'C')
                    # delta_n = cs_data[cs_list][residue]['N']-rc_shifts.get_value(residue[-1],'N')
                    # delta_h = cs_data[cs_list][residue]['H']-rc_shifts.get_value(residue[-1],'H')
                    delta_ca_cb = delta_ca-delta_cb
                    lacs_data[cs_list]['d_ca'].append(delta_ca)
                    lacs_data[cs_list]['d_cb'].append(delta_cb)
                    lacs_data[cs_list]['d_ca_cb'].append(delta_ca_cb)
                    lacs_data[cs_list]['tag'].append(residue)
                except KeyError:
                    pass
    fit_data ={}
    for cs_list in lacs_data:
        if cs_list not in fit_data:
            fit_data = {}
        x_p=[lacs_data[cs_list]['d_ca_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if lacs_data[cs_list]['d_ca_cb'][i] >= 0.0 ]
        y_p = [lacs_data[cs_list]['d_ca'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        x_n = [lacs_data[cs_list]['d_ca_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n = [lacs_data[cs_list]['d_ca'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        slop_p,offset_p = np.polyfit(x_p,y_p,1)
        slop_n,offset_n = np.polyfit(x_n,y_n,1)
        fit_data[cs_list] = {'slope_p':round(slop_p,2),'offset_p':round(offset_p,2),'slope_n':round(slop_n,2),'offset_n':round(offset_n,2),'p_max':max(x_p),'n_min':min(x_n)}


    fit_line={}
    for cs_list in fit_data:
        if cs_list not in fit_line:
            fit_line[cs_list]={}
        xp = np.linspace(0,fit_data[cs_list]['p_max'],100)
        xn = np.linspace(fit_data[cs_list]['n_min'],0,100)
        yp = fit_data[cs_list]['slope_p']*xp+fit_data[cs_list]['offset_p']
        yn = fit_data[cs_list]['slope_n']*xn+fit_data[cs_list]['offset_n']
        fit_line[cs_list] = {'xp':xp,'xn':xn,'yp':yp,'yn':yn}


    for cs_list in lacs_data:
        fig = px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_ca'],hover_name=lacs_data[cs_list]['tag'])
        fig.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_p']},Offset:{fit_data[cs_list]['offset_p']}',line=dict(color='red', dash='solid')))
        fig.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_n']},Offset:{fit_data[cs_list]['offset_n']}',line=dict(color='red', dash='dash')))
        fig.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$')
        fig.update_layout(yaxis_title=r'$\Delta\delta C^{\alpha}$')
        fig.show()

if __name__ == "__main__":
    lacs('../scratch/bmr52864_3.str')