#!/usr/bin/env python3
"""
This script is python version og LACS software  https://link.springer.com/article/10.1007/s10858-005-1717-0
"""
import pynmrstar
import plotly.express as px
import plotly.graph_objects as go
from random_coil import RandomCoil
import numpy as np
import statsmodels.api as sm

ONE_TO_THREE = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                               'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                               'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                               'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
THREE_TO_ONE = dict([(value, key) for key, value in ONE_TO_THREE.items()])


def read_str_file(file_name):
    """
    Reads a str file and returns a dictionary of chemical shifts
    :param file_name: NMR-STAR file name with path
    :return: Chemical shift dictionary
    """
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

def lacs(str_file,data_id = 'LACS_analysis',rc_model = None):
    """
    Calculate chemical shift offsets for a given NMR-STAR model with default random coil model
    :param str_file: NMR-STAR file
    :param rc_model: random coil chemical shift models
    :param rc_model: random coil chemical shift models
    :return: chemical shift offsets as dictionary
    """
    cs_data = read_str_file(str_file)
    rc_shifts = RandomCoil()
    atom_list = rc_shifts.atoms()
    lacs_data = {}
    for cs_list in cs_data:
        if cs_list not in lacs_data:
            lacs_data[cs_list] = {'d_ca':[],'d_cb':[],'d_ca_cb':[],'tag':[],'d_c':[],'d_n':[],'d_h':[],
                                  'd_c_x':[],'d_n_x':[],'d_h_x':[],'tag_c':[],'tag_n':[],'tag_h':[]}
        for residue in cs_data[cs_list]:
            if residue[-1] in THREE_TO_ONE:
                try:
                    delta_ca = cs_data[cs_list][residue]['CA']-rc_shifts.get_value(residue[-1],'CA',rc_name=rc_model)
                except KeyError:
                    print (f'Atom CA not found in {residue}')
                    delta_ca = None
                try:
                    delta_cb = cs_data[cs_list][residue]['CB']-rc_shifts.get_value(residue[-1],'CB',rc_name=rc_model)
                except KeyError:
                    print(f'Atom CB not found in {residue}')
                    delta_cb = None
                try:
                    delta_c = cs_data[cs_list][residue]['C']-rc_shifts.get_value(residue[-1],'C')
                except KeyError:
                    print(f'Atom C not found in {residue}')
                    delta_c = None
                try:
                    delta_n = cs_data[cs_list][residue]['N']-rc_shifts.get_value(residue[-1],'N')
                except KeyError:
                    print(f'Atom N not found in {residue}')
                    delta_n = None
                try:
                    delta_h = cs_data[cs_list][residue]['H']-rc_shifts.get_value(residue[-1],'H')
                except KeyError:
                    print(f'Atom H not found in {residue}')
                    delta_h = None
                if delta_ca is not None and delta_cb is not None:
                    delta_ca_cb = delta_ca-delta_cb
                    lacs_data[cs_list]['d_ca'].append(delta_ca)
                    lacs_data[cs_list]['d_cb'].append(delta_cb)
                    lacs_data[cs_list]['d_ca_cb'].append(delta_ca_cb)
                    lacs_data[cs_list]['tag'].append(residue)
                    if delta_c is not None:
                        lacs_data[cs_list]['d_c'].append(delta_c)
                        lacs_data[cs_list]['d_c_x'].append(delta_ca_cb)
                        lacs_data[cs_list]['tag_c'].append(residue)
                    if delta_n is not None:
                        lacs_data[cs_list]['d_n'].append(delta_n)
                        lacs_data[cs_list]['d_n_x'].append(delta_ca_cb)
                        lacs_data[cs_list]['tag_n'].append(residue)
                    if delta_h is not None:
                        lacs_data[cs_list]['d_h'].append(delta_h)
                        lacs_data[cs_list]['d_h_x'].append(delta_ca_cb)
                        lacs_data[cs_list]['tag_h'].append(residue)
    fit_data ={}
    for cs_list in lacs_data:
        if cs_list not in fit_data:
            fit_data[cs_list] = {}
        fit_data[cs_list]['ca'] = fit_data_rlm(lacs_data[cs_list]['d_ca_cb'], lacs_data[cs_list]['d_ca'],
                                           lacs_data[cs_list]['tag'])
        fit_data[cs_list]['cb'] = fit_data_rlm(lacs_data[cs_list]['d_ca_cb'], lacs_data[cs_list]['d_cb'],
                                           lacs_data[cs_list]['tag'])
        fit_data[cs_list]['c'] = fit_data_rlm(lacs_data[cs_list]['d_c_x'], lacs_data[cs_list]['d_c'],
                                           lacs_data[cs_list]['tag_c'])
        fit_data[cs_list]['n'] = fit_data_rlm(lacs_data[cs_list]['d_n_x'], lacs_data[cs_list]['d_n'],
                                           lacs_data[cs_list]['tag_n'])
        fit_data[cs_list]['h'] = fit_data_rlm(lacs_data[cs_list]['d_h_x'], lacs_data[cs_list]['d_h'],
                                           lacs_data[cs_list]['tag_h'])

    for cs_list in fit_data:
        for atom in fit_data[cs_list]:
            fig = px.scatter(x=fit_data[cs_list][atom]['x_n']+fit_data[cs_list][atom]['x_p'],
                             y = fit_data[cs_list][atom]['y_n']+fit_data[cs_list][atom]['y_p'],
                             hover_name=fit_data[cs_list][atom]['tag_n']+fit_data[cs_list][atom]['tag_p'])
            fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_p'],
                                     y=fit_data[cs_list][atom]['y_p'],
                                     mode='markers+text',
                                     text=fit_data[cs_list][atom]['tag_p'],
                                     textposition='top center',
                                     marker=dict(color=fit_data[cs_list][atom]['outliers_p'],size=10),showlegend=False))
            fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_n'],
                                     y=fit_data[cs_list][atom]['y_n'],
                                     mode='markers+text',
                                     text=fit_data[cs_list][atom]['tag_n'],
                                     textposition='top center',
                                     marker=dict(color=fit_data[cs_list][atom]['outliers_n'], size=10),
                                     showlegend=False))
            fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_p'],
                                        y=fit_data[cs_list][atom]['fittedvalues_p'],
                                        mode='lines',
                                        name=f'Slope:{fit_data[cs_list][atom]['slope_p']},Offset:{fit_data[cs_list][atom]['offset_p']} ',
                                        line=dict(color='green', dash='solid')))
            fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_n'],
                                     y=fit_data[cs_list][atom]['fittedvalues_n'],
                                     mode='lines',
                                     name=f'Slope:{fit_data[cs_list][atom]['slope_n']},Offset:{fit_data[cs_list][atom]['offset_n']} ',
                                     line=dict(color='green', dash='dash')))
            fig.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$')
            if atom == 'ca':
                fig.update_layout(yaxis_title=r'$\Delta\delta C^{\alpha}$')
            elif atom == 'cb':
                fig.update_layout(yaxis_title=r'$\Delta\delta C^{\beta}$')
            elif atom == 'c':
                fig.update_layout(yaxis_title=r'$\Delta\delta C$')
            elif atom == 'n':
                fig.update_layout(yaxis_title=r'$\Delta\delta N$')
            elif atom == 'h':
                fig.update_layout(yaxis_title=r'$\Delta\delta H $')
            else:
                raise ValueError(f'Atom {atom} not supported')
            fig.write_html(f'../scratch/{data_id}_{atom}.html',include_mathjax='cdn')
            fig.write_image(f'../scratch/{data_id}_{atom}.pdf')

            fig = px.bar(x=fit_data[cs_list][atom]['res_no'],
                         y=fit_data[cs_list][atom]['prob'],labels={'x':'Residue Number','y':'Probability of being outlier'})
            fig.update_layout(yaxis_range=[0, 1])
            fig.write_image(f'../scratch/{data_id}_{atom}_prob.pdf')
            fig.write_html(f'../scratch/{data_id}_{atom}_prob.html')
    offset={}
    for cs_list in fit_data:
        offset_ca = (fit_data[cs_list]['ca']['offset_p']+fit_data[cs_list]['ca']['offset_n'])/2.0
        offset_cb = (fit_data[cs_list]['cb']['offset_p'] + fit_data[cs_list]['cb']['offset_n']) / 2.0
        offset[cs_list]=(offset_ca,offset_cb)
    return offset



def fit_data_rlm(x,y,tag,outlier_cutoff = 0.9):
    x_p = [x[i] for i in range(len(x)) if x[i]>=0]
    y_p = [y[i] for i in range(len(x)) if x[i]>=0]
    x_n = [x[i] for i in range(len(x)) if x[i] < 0]
    y_n = [y[i] for i in range(len(x)) if x[i] < 0]
    tag_p = [f'{tag[i][-2]}{THREE_TO_ONE[tag[i][-1]]}' for i in range(len(x)) if x[i] >= 0.0]
    tag_n = [f'{tag[i][-2]}{THREE_TO_ONE[tag[i][-1]]}' for i in range(len(x)) if x[i] < 0.0]
    X_p = sm.add_constant(x_p)
    X_n = sm.add_constant(x_n)
    rlm_p = sm.RLM(y_p, X_p, M=sm.robust.norms.TukeyBiweight())
    rlm_n = sm.RLM(y_n, X_n, M=sm.robust.norms.TukeyBiweight())
    rlm_p_result = rlm_p.fit()
    rlm_n_result = rlm_n.fit()
    offset_p,slope_p = round(rlm_p_result.params[0],2),round(rlm_p_result.params[1],2)
    offset_n, slope_n = round(rlm_n_result.params[0],2), round(rlm_n_result.params[1],2)
    prob_p = [1.0-i for i in rlm_p_result.weights.tolist()]
    prob_n = [1.0-i for i in rlm_n_result.weights.tolist()]
    outliers_p = ['green' if i > outlier_cutoff else 'red' for i in prob_p]
    outliers_n = ['green' if i > outlier_cutoff else 'red' for i in prob_n]
    fittedvalues_p = rlm_p_result.fittedvalues.tolist()
    fittedvalues_n = rlm_n_result.fittedvalues.tolist()
    res_no=[]
    prob=[]
    for i in range(len(tag_p)):
        res_no.append(int(tag_p[i][:-1]))
        prob.append(prob_p[i])
    for i in range(len(tag_n)):
        res_no.append(int(tag_n[i][:-1]))
        prob.append(prob_n[i])
    fit_data = {'x_p':x_p,
                'x_n':x_n,
                'y_p':y_p,
                'y_n':y_n,
                'tag_p':tag_p,
                'tag_n':tag_n,
                'slope_p':slope_p,
                'slope_n':slope_n,
                'offset_p':offset_p,
                'offset_n':offset_n,
                'outliers_p':outliers_p,
                'outliers_n':outliers_n,
                'fittedvalues_p':fittedvalues_p,
                'fittedvalues_n':fittedvalues_n,
                'res_no':res_no,
                'prob':prob}
    return fit_data

if __name__ == "__main__":
    print (lacs('../scratch/bmr25421_3.str',data_id='25421'))