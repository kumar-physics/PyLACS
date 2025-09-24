import csv
from random_coil import RandomCoil
import plotly.express as px
from lacs import run_lacs2
import pandas as pd

def read_csv(csvfile):
    rc = RandomCoil()
    with open(csvfile, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    json_data ={}
    for row in data:
        if row[-1] not in json_data:
            json_data[row[-1]]={}
        if row[3] not in json_data[row[-1]]:
            json_data[row[-1]][row[3]] = {}
        if row[4] not in json_data[row[-1]][row[3]]:
            json_data[row[-1]][row[3]][row[4]] ={}
        if int(row[0]) not in json_data[row[-1]][row[3]][row[4]]:
            json_data[row[-1]][row[3]][row[4]][int(row[0])]= {}
        try:
            json_data[row[-1]][row[3]][row[4]][int(row[0])][row[2]]= (float(row[6]),round(rc.get_value(row[1],row[2])-float(row[6]),3),row[7],row[1])
        except ValueError:
            pass
    return json_data

def read_csv2(csvfile):
    rc = RandomCoil()
    with open(csvfile, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    json_data ={}
    for row in data:
        key = (row[-1],row[3],row[4],row[0],row[1])
        if key not in json_data:
            json_data[key]={}
        if row[2] in ['CA','CB','C','HA','N']:
            try:
                json_data[key][row[2]]= float(row[6])
            except ValueError:
                pass
    return {'1':json_data}

def plot_data(csvfile):
    data = read_csv(csvfile)
    y=[]
    x=[]
    ss=[]
    tag=[]
    r=[]
    for ent in data:
        for lst in data[ent]:
            for ch in data[ent][lst]:
                for res in sorted(data[ent][lst][ch]):
                    try:
                        if abs(data[ent][lst][ch][res]['CA'][1])<50 and abs(data[ent][lst][ch][res]['CB'][1])<50:# and abs(data[ent][lst][ch][res]['HA'][1])<10:
                            d= data[ent][lst][ch][res]['CA'][1]-data[ent][lst][ch][res]['CB'][1]
                            s= data[ent][lst][ch][res]['CA'][1]#+data[ent][lst][ch][res+1]['CA'][1]
                            x.append(d)
                            y.append(s)
                            r.append(data[ent][lst][ch][res]['CA'][3])
                            ss.append(data[ent][lst][ch][res]['CA'][2])
                            tag.append(ent)
                    except KeyError:
                        pass


    fig = px.scatter(x=x, y=y,color=ss,facet_col=r,facet_col_wrap=4,hover_name=tag,labels={'x':"ΔδCA − ΔδCB",'y':"ΔδC4A"})#, marginal_x='histogram', marginal_y='histogram',labels={'x':"ΔδCA − ΔδCB",'y':"ΔδHA"})
    fig.write_html('CA_SS3.html')

def lacs(fname):
    cs = read_csv2(fname)
    run_lacs2(cs)
if __name__ == "__main__":
    #plot_data('../scratch/psb_bmrb.csv')
    lacs('../scratch/psb_bmrb.csv')