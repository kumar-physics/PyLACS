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
import csv
import pandas as pd

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
    try:
        ent = pynmrstar.Entry.from_file(file_name)
    except FileNotFoundError:
        return {}
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
    if len(cs_data) == 0:
       return {}
    rc_shifts = RandomCoil()
    atom_list = rc_shifts.atoms()
    lacs_data = {}
    for cs_list in cs_data:
        if cs_list not in lacs_data:
            lacs_data[cs_list] = {'d_ca':[],'d_cb':[],'d_ca_cb':[],'d_c':[],'d_n':[],'d_h':[],
                                  'd_c_x':[],'d_n_x':[],'d_h_x':[],'tag_c':[],'tag_n':[],'tag_h':[],'tag_ca':[],'tag_cb':[]}
        for residue in cs_data[cs_list]:
            if residue[-1] in THREE_TO_ONE:
                try:
                    delta_ca = cs_data[cs_list][residue]['CA']-rc_shifts.get_value(residue[-1],'CA',rc_name=rc_model)
                except KeyError:
                    # print (f'Atom CA not found in {residue}')
                    delta_ca = None
                try:
                    delta_cb = cs_data[cs_list][residue]['CB']-rc_shifts.get_value(residue[-1],'CB',rc_name=rc_model)
                except KeyError:
                    # print(f'Atom CB not found in {residue}')
                    delta_cb = None
                try:
                    delta_c = cs_data[cs_list][residue]['C']-rc_shifts.get_value(residue[-1],'C',rc_name=rc_model)
                except KeyError:
                    # print(f'Atom C not found in {residue}')
                    delta_c = None
                try:
                    delta_n = cs_data[cs_list][residue]['N']-rc_shifts.get_value(residue[-1],'N',rc_name=rc_model)
                except KeyError:
                    # print(f'Atom N not found in {residue}')
                    delta_n = None
                try:
                    delta_h = cs_data[cs_list][residue]['H']-rc_shifts.get_value(residue[-1],'H',rc_name=rc_model)
                except KeyError:
                    # print(f'Atom H not found in {residue}')
                    delta_h = None
                if delta_ca is not None and delta_cb is not None:
                    delta_ca_cb = delta_ca-delta_cb
                    lacs_data[cs_list]['d_ca'].append(delta_ca)
                    lacs_data[cs_list]['d_cb'].append(delta_cb)
                    lacs_data[cs_list]['d_ca_cb'].append(delta_ca_cb)
                    lacs_data[cs_list]['tag_ca'].append(residue)
                    lacs_data[cs_list]['tag_cb'].append(residue)
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
    atoms = ['ca','cb']
    fit_data ={}
    for cs_list in lacs_data:
        if len(lacs_data[cs_list]['d_c']):atoms.append('c')
        if len(lacs_data[cs_list]['d_n']): atoms.append('n')
        if len(lacs_data[cs_list]['d_h']): atoms.append('h')
        if cs_list not in fit_data:
            fit_data[cs_list] = {}
        for atom in atoms:
            if atom not in ['ca','cb']:
                fit_data[cs_list][atom] = fit_data_rlm(lacs_data[cs_list][f'd_{atom}_x'], lacs_data[cs_list][f'd_{atom}'],
                                                       lacs_data[cs_list][f'tag_{atom}'])
            else:
                fit_data[cs_list][atom] = fit_data_rlm(lacs_data[cs_list]['d_ca_cb'], lacs_data[cs_list][f'd_{atom}'],
                                           lacs_data[cs_list][f'tag_{atom}'])


    # for cs_list in fit_data:
    #     for atom in fit_data[cs_list]:
    #         fig = px.scatter(x=fit_data[cs_list][atom]['x_n']+fit_data[cs_list][atom]['x_p'],
    #                          y = fit_data[cs_list][atom]['y_n']+fit_data[cs_list][atom]['y_p'],
    #                          hover_name=fit_data[cs_list][atom]['tag_n']+fit_data[cs_list][atom]['tag_p'])
    #         fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_p'],
    #                                  y=fit_data[cs_list][atom]['y_p'],
    #                                  mode='markers+text',
    #                                  text=fit_data[cs_list][atom]['tag_p'],
    #                                  textposition='top center',
    #                                  marker=dict(color=fit_data[cs_list][atom]['outliers_p'],size=10),showlegend=False))
    #         fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_n'],
    #                                  y=fit_data[cs_list][atom]['y_n'],
    #                                  mode='markers+text',
    #                                  text=fit_data[cs_list][atom]['tag_n'],
    #                                  textposition='top center',
    #                                  marker=dict(color=fit_data[cs_list][atom]['outliers_n'], size=10),
    #                                  showlegend=False))
    #         fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_p'],
    #                                     y=fit_data[cs_list][atom]['fittedvalues_p'],
    #                                     mode='lines',
    #                                     name=f'Slope:{fit_data[cs_list][atom]['slope_p']},Offset:{fit_data[cs_list][atom]['offset_p']} ',
    #                                     line=dict(color='green', dash='solid')))
    #         fig.add_trace(go.Scatter(x=fit_data[cs_list][atom]['x_n'],
    #                                  y=fit_data[cs_list][atom]['fittedvalues_n'],
    #                                  mode='lines',
    #                                  name=f'Slope:{fit_data[cs_list][atom]['slope_n']},Offset:{fit_data[cs_list][atom]['offset_n']} ',
    #                                  line=dict(color='green', dash='dash')))
    #         fig.update_layout(xaxis_title=r'$\Delta\delta C^{\alpha}-\Delta\delta C^{\beta}$')
    #         if atom == 'ca':
    #             fig.update_layout(yaxis_title=r'$\Delta\delta C^{\alpha}$')
    #         elif atom == 'cb':
    #             fig.update_layout(yaxis_title=r'$\Delta\delta C^{\beta}$')
    #         elif atom == 'c':
    #             fig.update_layout(yaxis_title=r'$\Delta\delta C$')
    #         elif atom == 'n':
    #             fig.update_layout(yaxis_title=r'$\Delta\delta N$')
    #         elif atom == 'h':
    #             fig.update_layout(yaxis_title=r'$\Delta\delta H $')
    #         else:
    #             raise ValueError(f'Atom {atom} not supported')
    #         fig.write_html(f'../scratch/{data_id}_{atom}.html',include_mathjax='cdn')
    #         fig.write_image(f'../scratch/{data_id}_{atom}.pdf')
    #
    #         fig = px.bar(x=fit_data[cs_list][atom]['res_no'],
    #                      y=fit_data[cs_list][atom]['prob'],labels={'x':'Residue Number','y':'Probability of being outlier'})
    #         fig.update_layout(yaxis_range=[0, 1])
    #         fig.write_image(f'../scratch/{data_id}_{atom}_prob.pdf')
    #         fig.write_html(f'../scratch/{data_id}_{atom}_prob.html')
    offset={}
    for cs_list in fit_data:
        if cs_list not in offset:
            offset[cs_list]={}
        for atom in fit_data[cs_list]:
            offset[cs_list][atom]=round((fit_data[cs_list][atom]['offset_p']+fit_data[cs_list][atom]['offset_n'])/2.0,2)
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

def read_csv(csv_file):
    flg=False
    bid_ca=[]
    bid_c=[]
    bid_h=[]
    tag_ca=[]
    tag_c=[]
    tag_h=[]
    CA=[]
    C=[]
    H=[]
    f=open('Pylacs_results.csv','w')
    with open (csv_file,'r') as csvfile:
        csvFile = csv.reader(csvfile)
        for row in csvFile:
            bmrb_id = row[0]
            if bmrb_id not in ['4094','4109','4141','4146','4165','4402','4813','4898',
                                       '4955','4994','5009','5226','5349','5361','5566','5571',
                                       '5623','5753']:
                bmrb_file = f'/reboxitory/2025/06/BMRB/macromolecules/bmr{bmrb_id}/bmr{bmrb_id}_3.str'
                offset = lacs(bmrb_file,bmrb_id)
                if len(offset)==9:
                    ca=''
                    c=''
                    h=''
                else:
                    try:
                        ca = offset['1']['ca']
                    except KeyError:
                        ca = ''
                    try:
                        c = offset['1']['c']
                    except KeyError:
                        c=''
                    try:
                        h = offset['1']['h']
                    except KeyError:
                        h=''


                    if ca!='':
                        bid_ca.append(bmrb_id)
                        CA.append(ca)
                        tag_ca.append('PyLACS')
                        bid_ca.append(bmrb_id)
                        if row[4]=='':
                            CA.append(None)
                        else:
                            CA.append(float(row[4]))
                        tag_ca.append('LACS')
                        bid_ca.append(bmrb_id)
                        if row[4] == '':
                            CA.append(None)
                        else:
                            CA.append(float(row[5]))
                        tag_ca.append('RefBD')
                    if c!='':
                        bid_c.append(bmrb_id)
                        C.append(ca)
                        tag_c.append('PyLACS')
                        bid_c.append(bmrb_id)
                        if row[6]=='':
                            C.append(None)
                        else:
                            C.append(float(row[6]))
                        tag_c.append('LACS')
                        bid_c.append(bmrb_id)
                        if row[7] == '':
                            C.append(None)
                        else:
                            C.append(float(row[7]))
                        tag_c.append('RefBD')
                    if h!='':
                        bid_h.append(bmrb_id)
                        H.append(ca)
                        tag_h.append('PyLACS')
                        bid_h.append(bmrb_id)
                        if row[8]=='':
                            H.append(None)
                        else:
                            H.append(float(row[8]))
                        tag_h.append('LACS')
                        bid_h.append(bmrb_id)
                        if row[9] == '':
                            H.append(None)
                        else:
                            H.append(float(row[9]))
                        tag_h.append('RefBD')
                    f.write(f'{row[0]},{row[1]},{row[2]},{row[3]},{ca},{row[4]},{row[5]},{c},{row[6]},{row[7]},{h},{row[8]},{row[9]}\n')

        fig.show()

def plot_data(csv_file):
    data=[]
    with open (csv_file,'r') as csvfile:
        csvFile = csv.reader(csvfile)
        data = [list(row) for row in zip(*list(csvFile))]
        bmrb_id = data[0]
        pylacs_ca =[-float(i) if i !='' else None for i in data[4]]
        lacs_ca = [float(i) if i !='' else None for i in  data[5]]
        refdb_ca = [float(i) if i !='' else None for i in data[6]]
    df = pd.DataFrame({
        'BMRB_ID':bmrb_id,
        'PyLACS':pylacs_ca,
        'LACS':lacs_ca,
        'RefDB':refdb_ca
    })
    df2 = pd.DataFrame({
        'BMRB_ID': bmrb_id,
        'RefDB':refdb_ca,
        'PyLACS': [refdb_ca[i]-pylacs_ca[i] if pylacs_ca[i] is not None else None for i in range(len(refdb_ca))],
        'LACS': [refdb_ca[i]-lacs_ca[i] if lacs_ca[i] is not None else None for i in range(len(refdb_ca))]
    })

    fig = px.scatter(df,x='RefDB',y=['LACS','PyLACS'],hover_name='BMRB_ID')
    fig.write_html('RefDB_vs_LACS.html')
    fig=px.bar(df2,x='BMRB_ID',y=['PyLACS','LACS'],barmode='group')
    fig.write_html('RefDB_vs_LACS_diff_bar.html')
    fig = px.scatter(df2,x='PyLACS',y='LACS',hover_name='BMRB_ID')
    fig.write_html('RefDB_vs_LACS_correlation.html')
    fig = px.scatter(df2,x='RefDB',y=['PyLACS','LACS'],hover_name='BMRB_ID')
    fig.write_html('RefDB_vs_LACS_correlation2.html')




if __name__ == "__main__":
    plot_data('Pylacs_results.csv')
    #print (lacs('../scratch/bmr4998_3.str',data_id=4998))
    #read_csv('Table_S1_clean.csv')
    # bmrb_id = 4109
    # file_name = f'/reboxitory/2025/06/BMRB/macromolecules/bmr{bmrb_id}/bmr{bmrb_id}_3.str'
    # print(lacs(str_file=file_name,data_id=bmrb_id))