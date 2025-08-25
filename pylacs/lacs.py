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
            lacs_data[cs_list] = {'d_ca':[],'d_cb':[],'d_ca_cb':[],'tag':[],'d_c':[],'d_n':[],'d_h':[]}
        for residue in cs_data[cs_list]:
            if residue[-1] in THREE_TO_ONE:
                try:
                    delta_ca = cs_data[cs_list][residue]['CA']-rc_shifts.get_value(residue[-1],'CA',rc_name=rc_model)
                except KeyError:
                    print (f'Atom CA not found in {residue}')
                    delta_ca = np.nan
                try:
                    delta_cb = cs_data[cs_list][residue]['CB']-rc_shifts.get_value(residue[-1],'CB',rc_name=rc_model)
                except KeyError:
                    print(f'Atom CB not found in {residue}')
                    delta_cb = np.nan
                try:
                    delta_c = cs_data[cs_list][residue]['C']-rc_shifts.get_value(residue[-1],'C')
                except KeyError:
                    print(f'Atom C not found in {residue}')
                    delta_c = np.nan
                try:
                    delta_n = cs_data[cs_list][residue]['N']-rc_shifts.get_value(residue[-1],'N')
                except KeyError:
                    print(f'Atom N not found in {residue}')
                    delta_n = np.nan
                try:
                    delta_h = cs_data[cs_list][residue]['H']-rc_shifts.get_value(residue[-1],'H')
                except KeyError:
                    print(f'Atom H not found in {residue}')
                    delta_h = np.nan
                if delta_ca is not None and delta_cb is not None:
                    delta_ca_cb = delta_ca-delta_cb
                    lacs_data[cs_list]['d_ca'].append(delta_ca)
                    lacs_data[cs_list]['d_cb'].append(delta_cb)
                    lacs_data[cs_list]['d_ca_cb'].append(delta_ca_cb)
                    lacs_data[cs_list]['d_c'].append(delta_c)
                    lacs_data[cs_list]['d_n'].append(delta_n)
                    lacs_data[cs_list]['d_h'].append(delta_h)
                    lacs_data[cs_list]['tag'].append(residue)
    fit_data ={}
    for cs_list in lacs_data:
        if cs_list not in fit_data:
            fit_data = {}
        x_p_ca_cb=[lacs_data[cs_list]['d_ca_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
                   lacs_data[cs_list]['d_ca_cb'][i] >= 0.0 ]
        y_p_ca = [lacs_data[cs_list]['d_ca'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        y_p_cb = [lacs_data[cs_list]['d_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        y_p_c = [lacs_data[cs_list]['d_c'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        y_p_n = [lacs_data[cs_list]['d_n'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        y_p_h = [lacs_data[cs_list]['d_h'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] >= 0.0]
        x_n_ca_cb = [lacs_data[cs_list]['d_ca_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n_ca = [lacs_data[cs_list]['d_ca'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n_cb = [lacs_data[cs_list]['d_cb'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n_c = [lacs_data[cs_list]['d_c'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n_n = [lacs_data[cs_list]['d_n'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        y_n_h = [lacs_data[cs_list]['d_h'][i] for i in range(len(lacs_data[cs_list]['d_ca_cb'])) if
               lacs_data[cs_list]['d_ca_cb'][i] < 0.0]
        slop_p_ca,offset_p_ca = np.polyfit(x_p_ca_cb,y_p_ca,1)
        slop_n_ca,offset_n_ca = np.polyfit(x_n_ca_cb,y_n_ca,1)
        slop_p_cb, offset_p_cb = np.polyfit(x_p_ca_cb, y_p_cb, 1)
        slop_n_cb, offset_n_cb = np.polyfit(x_n_ca_cb, y_n_cb, 1)
        slop_p_c, offset_p_c = np.polyfit(x_p_ca_cb, y_p_c, 1)
        slop_n_c, offset_n_c = np.polyfit(x_n_ca_cb, y_n_c, 1)
        slop_p_n, offset_p_n = np.polyfit(x_p_ca_cb, y_p_n, 1)
        slop_n_n, offset_n_n = np.polyfit(x_n_ca_cb, y_n_n, 1)
        slop_p_h, offset_p_h = np.polyfit(x_p_ca_cb, y_p_h, 1)
        slop_n_h, offset_n_h = np.polyfit(x_n_ca_cb, y_n_h, 1)
        fit_data[cs_list] = {'slope_p_ca':round(slop_p_ca,2),'offset_p_ca':round(offset_p_ca,2),'slope_n_ca':round(slop_n_ca,2),'offset_n_ca':round(offset_n_ca,2),
                             'slope_p_cb':round(slop_p_cb,2),'offset_p_cb':round(offset_p_cb,2),'slope_n_cb':round(slop_n_cb,2),'offset_n_cb':round(offset_n_cb,2),
                             'slope_p_c':round(slop_p_c,2),'offset_p_c':round(offset_p_c,2),'slope_n_c':round(slop_n_c,2),'offset_n_c':round(offset_n_c,2),
                             'slope_p_n':round(slop_p_n,2),'offset_p_n':round(offset_p_n,2),'slope_n_n':round(slop_n_n,2),'offset_n_n':round(offset_n_n,2),
                             'slope_p_h':round(slop_p_h,2),'offset_p_h':round(offset_p_h,2),'slope_n_h':round(slop_n_h,2),'offset_n_h':round(offset_n_h,2),
                             'p_max':max(x_p_ca_cb),'n_min':min(x_n_ca_cb)}


    fit_line={}
    for cs_list in fit_data:
        if cs_list not in fit_line:
            fit_line[cs_list]={}
        xp = np.linspace(0,fit_data[cs_list]['p_max'],100)
        xn = np.linspace(fit_data[cs_list]['n_min'],0,100)
        yp_ca = fit_data[cs_list]['slope_p_ca']*xp+fit_data[cs_list]['offset_p_ca']
        yn_ca = fit_data[cs_list]['slope_n_ca']*xn+fit_data[cs_list]['offset_n_ca']
        yp_cb = fit_data[cs_list]['slope_p_cb']*xp+fit_data[cs_list]['offset_p_cb']
        yn_cb = fit_data[cs_list]['slope_n_cb']*xn+fit_data[cs_list]['offset_n_cb']
        yp_c = fit_data[cs_list]['slope_p_c']*xp+fit_data[cs_list]['offset_p_c']
        yn_c = fit_data[cs_list]['slope_n_c']*xn+fit_data[cs_list]['offset_n_c']
        yp_n = fit_data[cs_list]['slope_p_n']*xp+fit_data[cs_list]['offset_p_n']
        yn_n = fit_data[cs_list]['slope_n_n']*xn+fit_data[cs_list]['offset_n_n']
        yp_h = fit_data[cs_list]['slope_p_h']*xp+fit_data[cs_list]['offset_p_h']
        yn_h = fit_data[cs_list]['slope_n_h']*xn+fit_data[cs_list]['offset_n_h']
        fit_line[cs_list] = {'xp':xp,'xn':xn,
                             'yp_ca':yp_ca,'yn_ca':yn_ca,
                             'yp_cb':yp_cb,'yn_cb':yn_cb,
                             'yp_c':yp_c,'yn_c':yn_c,
                             'yp_n':yp_n,'yn_n':yn_n,
                             'yp_h':yp_h,'yn_h':yn_h}


    for cs_list in lacs_data:
        fig_ca = px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_ca'],hover_name=lacs_data[cs_list]['tag'])
        fig_ca.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp_ca'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_p_ca']},Offset:{fit_data[cs_list]['offset_p_ca']}',line=dict(color='red', dash='solid')))
        fig_ca.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn_ca'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_n_ca']},Offset:{fit_data[cs_list]['offset_n_ca']}',line=dict(color='red', dash='dash')))
        fig_ca.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$')
        fig_ca.update_layout(yaxis_title=r'$\Delta\delta C^{\alpha}$')
        fig_ca.show()
        fig_cb = px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_cb'],hover_name=lacs_data[cs_list]['tag'])
        fig_cb.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp_cb'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_p_cb']},Offset:{fit_data[cs_list]['offset_p_cb']}',line=dict(color='red', dash='solid')))
        fig_cb.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn_cb'],mode='lines',name=f'Slope:{fit_data[cs_list]['slope_n_cb']},Offset:{fit_data[cs_list]['offset_n_cb']}',line=dict(color='red', dash='dash')))
        fig_cb.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$')
        fig_cb.update_layout(yaxis_title=r'$\Delta\delta C^{\beta}$')
        fig_cb.show()
        fig_c=px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_c'],hover_name=lacs_data[cs_list]['tag'])
        fig_c.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp_c'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_p_c']}, <b>Offset:</b> {fit_data[cs_list]['offset_p_c']}',line=dict(color='red', dash='solid')))
        fig_c.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn_c'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_n_c']}, <b>Offset:</b> {fit_data[cs_list]['offset_n_c']}',line=dict(color='red', dash='dash')))
        fig_c.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$',yaxis_title=r'$\Delta\delta C$')
        fig_c.show()
        fig_n=px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_n'],hover_name=lacs_data[cs_list]['tag'])
        fig_n.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp_n'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_p_n']}, <b>Offset:</b> {fit_data[cs_list]['offset_p_n']}',line=dict(color='red', dash='solid')))
        fig_n.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn_n'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_n_n']}, <b>Offset:</b> {fit_data[cs_list]['offset_n_n']}',line=dict(color='red', dash='dash')))
        fig_n.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$',yaxis_title=r'$\Delta\delta N$')
        fig_n.show()
        fig_h=px.scatter(x=lacs_data[cs_list]['d_ca_cb'], y=lacs_data[cs_list]['d_h'],hover_name=lacs_data[cs_list]['tag'])
        fig_h.add_trace(go.Scatter(x=fit_line[cs_list]['xp'],y=fit_line[cs_list]['yp_h'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_p_h']}, <b>Offset:</b> {fit_data[cs_list]['offset_p_h']}',line=dict(color='red', dash='solid')))
        fig_h.add_trace(go.Scatter(x=fit_line[cs_list]['xn'],y=fit_line[cs_list]['yn_h'],mode='lines',name=f'<b>Slope:</b> {fit_data[cs_list]['slope_n_h']}, <b>Offset:</b> {fit_data[cs_list]['offset_n_h']}',line=dict(color='red', dash='dash')))
        fig_h.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$',yaxis_title=r'$\Delta\delta H$')
        fig_h.show()

if __name__ == "__main__":
    lacs('../scratch/bmr30196_3.str')