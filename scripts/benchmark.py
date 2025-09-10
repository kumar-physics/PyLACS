import json
import csv
import pandas as pd
import plotly.express as px
import numpy as np


def read_csv(panav_file):
    panav_err ={}
    with open(panav_file,'r') as panvFile:
        csvReader = csv.reader(panvFile)
        for row in csvReader:
            if int(row[-1])>24 and row[0] not in panav_err:
                panav_err[row[0]]={}
            if int(row[-1])>24:
                panav_err[row[0]][row[2]]=float(row[3])
    return panav_err

def read_json(lacs_file):
    with open(lacs_file, "r", encoding="utf-8") as f:
        lacs_err = json.load(f)
    return lacs_err

def compar_offsets(panav_err,lacs_err):
    bmrb_id =[]
    panav_ca=[]
    panav_cb=[]
    panav_c=[]
    panav_n=[]
    tukey_ca=[]
    tukey_cb=[]
    tukey_c=[]
    tukey_n=[]
    ran_ca=[]
    ran_cb=[]
    ran_c=[]
    ran_n=[]
    q_ca=[]
    q_cb=[]
    q_c=[]
    q_n=[]
    t_ca=[]
    t_cb=[]
    t_c=[]
    t_n=[]
    b_ca=[]
    b_cb=[]
    b_c=[]
    b_n=[]
    med_ca=[]

    panav = read_csv(panav_err)
    lacs = read_json(lacs_err)
    for entry_id in panav:
        ca=[]
        try:
            print (lacs[entry_id]['1'])
            bmrb_id.append(entry_id)
            try:
                panav_ca.append(panav[entry_id]['CA'])
            except KeyError:
                panav_ca.append(None)
            try:
                # panav_ca.append(panav[entry_id]['CB'])
                # panav_ca.append(panav[entry_id]['C'])
                # panav_ca.append(panav[entry_id]['N'])
                tukey_ca.append(lacs[entry_id]['1']['tukey']['ca'])
                if lacs[entry_id]['1']['tukey']['ca'] is not None: ca.append(lacs[entry_id]['1']['tukey']['ca'])
            except KeyError:
                tukey_ca.append(None)

            try:
                tukey_cb.append(lacs[entry_id]['1']['tukey']['cb'])
            except KeyError:
                tukey_cb.append(None)
            try:
                tukey_c.append(lacs[entry_id]['1']['tukey']['c'])
            except KeyError:
                tukey_c.append(None)
            try:
                tukey_n.append(lacs[entry_id]['1']['tukey']['n'])
            except KeyError:
                tukey_n.append(None)
            try:
                ran_ca.append(lacs[entry_id]['1']['ransac']['ca'])
                if lacs[entry_id]['1']['ransac']['ca'] is not None: ca.append(lacs[entry_id]['1']['ransac']['ca'])
            except KeyError:
                ran_ca.append(None)
            try:
                ran_cb.append(lacs[entry_id]['1']['ransac']['cb'])
            except KeyError:
                ran_cb.append(None)
            try:
                ran_c.append(lacs[entry_id]['1']['ransac']['c'])
            except KeyError:
                ran_c.append(None)
            try:
                ran_n.append(lacs[entry_id]['1']['ransac']['n'])
            except KeyError:
                ran_n.append(None)
            try:
                q_ca.append(lacs[entry_id]['1']['quantile']['ca'])
                if lacs[entry_id]['1']['quantile']['ca'] is not None: ca.append(lacs[entry_id]['1']['quantile']['ca'])
            except KeyError:
                q_ca.append(None)
            try:
                q_cb.append(lacs[entry_id]['1']['quantile']['cb'])
            except KeyError:
                q_cb.append(None)
            try:
                q_c.append(lacs[entry_id]['1']['quantile']['c'])
            except KeyError:
                q_c.append(None)
            try:
                q_n.append(lacs[entry_id]['1']['quantile']['n'])
            except KeyError:
                q_n.append(None)
            try:
                t_ca.append(lacs[entry_id]['1']['theilsen']['ca'])
                if lacs[entry_id]['1']['theilsen']['ca'] is not None: ca.append(lacs[entry_id]['1']['theilsen']['ca'])
            except KeyError:
                t_ca.append(None)
            try:
                t_cb.append(lacs[entry_id]['1']['theilsen']['cb'])
            except KeyError:
                t_cb.append(None)
            try:
                t_c.append(lacs[entry_id]['1']['theilsen']['c'])
            except KeyError:
                t_c.append(None)
            try:
                t_n.append(lacs[entry_id]['1']['bayes']['n'])
            except KeyError:
                t_n.append(None)
            try:
                b_ca.append(lacs[entry_id]['1']['bayes']['ca'])
                if lacs[entry_id]['1']['bayes']['ca'] is not None: ca.append(lacs[entry_id]['1']['bayes']['ca'])
            except KeyError:
                b_ca.append(None)
            try:
                b_cb.append(lacs[entry_id]['1']['bayes']['cb'])
            except KeyError:
                b_cb.append(None)
            try:
                b_c.append(lacs[entry_id]['1']['bayes']['c'])
            except KeyError:
                b_c.append(None)
            try:
                b_n.append(lacs[entry_id]['1']['bayes']['n'])
            except KeyError:
                b_n.append(None)
            print ("here",ca)
            if len(ca):
                med_ca.append(np.median(ca))
            else:
                med_ca.append(None)

        except KeyError:
            print (entry_id)
    print (len(bmrb_id),len(panav_ca),len(tukey_ca),len(t_ca),len(ran_ca),len(q_ca),len(b_ca))
    df = pd.DataFrame({  # ["tukey", "theilsen", "ransac", "quantile", "bayes"]
        'BMRB_ID': bmrb_id,
        'LACS(CA)':med_ca,
        'PANAV(CA)': panav_ca,
        # 'PANAV(C)': panav_c,
        # 'PANAV(CB)': panav_cb,
        # 'PANAV(N)' : panav_n,
        'Tukey(CA)': tukey_ca,
        # 'Tukey(C)': tukey_c,
        # 'Tukey(CB)': tukey_cb,
        # 'Tukey(N)' : tukey_n,
        'Theilsen(CA)': t_ca,
        # 'Theilsen(C)': t_c,
        # 'Theilsen(CB)': t_cb,
        # 'Theilsen(N)': t_n,
        'Ransac(CA)': ran_ca,
        # 'Ransac(C)': ran_c,
        # 'Ransac(CB)': ran_cb,
        # 'Ransac(N)': ran_n,
        'Quantile(CA)': q_ca,
        # 'Quantile(C)': q_c,
        # 'Quantile(CB)': q_cb,
        # 'Quantile(N)': q_n,
        'Bayes(CA)': b_ca,
        # 'Bayes(C)': b_c,
        # 'Bayes(CB)': b_cb,
        # 'Bayes(N)': b_n,


    })
    fig = px.scatter(df, x='PANAV(CA)',
                     y='LACS(CA)',
                     hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-30,  # Starting x-coordinate of the line
        y0=-30,  # Starting y-coordinate of the line
        x1=30,  # Ending x-coordinate of the line
        y1=30,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.show()
    fig.write_html('benchmark.html')
if __name__ == "__main__":
   compar_offsets('pannav_errors.csv','bmrb_offsets.json')