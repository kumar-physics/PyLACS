import csv
import argparse
import json
import math
from pathlib import Path
from statistics import mean, variance, median
from typing import Any, Dict, Iterable, Tuple
import numpy as np
import pandas as pd
import plotly.express as px
from matplotlib.pyplot import ylabel

FITS = ['bayes','theilsen','ransac','quantile','tukey']
RCS = ['avg','luk','pou','sch','wan','wis']
LACS_PATH = '/Users/kumaranbaskaran/lacs/'
#LACS_PATH = '/Users/kumaranbaskaran/lacs_refdb_Feb28/'
#FIGS_PATH = '/Users/kumaranbaskaran/Documents/PyLACS_Manuscript/figures/LACS_tableS1_283_entreis'
FIGS_PATH = '/Users/kumaranbaskaran/Documents/PyLACS_Manuscript/figures/RefDB_2426_entries'
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
    out=[]
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
                lacs_data = read_json_dict(lacs_output,fit)
                if lacs_data is not None:
                    for atm in lacs_data:
                        if atm == 'c':
                            atm2='CO'
                        else:
                            atm2=atm
                        if fit == 'bayes':
                            out.append([bmrb_id,atm,fit,rc,row[data[0].index(f'{atm2.upper()}_CSDIFF')],lacs_data[atm]['mean'],lacs_data[atm]['sd']])
                        else:
                            out.append([bmrb_id,atm, fit,rc, row[data[0].index(f'{atm2.upper()}_CSDIFF')],lacs_data[atm],0.0])
    with open('RefDB_Lacs_benchmark.csv', 'w') as f:
        csv_write =csv.writer(f,delimiter=',')
        csv_write.writerows(out)

def plot_data(csv_file):
    with open(csv_file) as f:
        data = []
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            data.append(row)

    for i in range(len(data)):
        try:
            data[i][4]=float(data[i][4])
        except ValueError:
            data[i][4] = float('nan')
        try:
            data[i][5] = float(data[i][5])
        except ValueError:
            data[i][5] = float('nan')
        try:
            data[i][6] = float(data[i][6])
        except ValueError:
            data[i][6] = float('nan')
        try:
            data[i].append(data[i][4]-data[i][5])
            data[i].append(abs(data[i][4] - data[i][5]))
        except ValueError:
            data[i].append(float('nan'))
            data[i].append(float('nan'))
    # data2= [i for i in data if i[1] == 'ca']
    # print (data[0])
    err_data=[]
    of=[]
    plot_atom=None
    for atm in ['n']:#,'cb','c','n','ha']:
        plot_atom=atm.upper()
        for rc in RCS:
            for fit in FITS:
                err=[]
                err_abs=[]
                for i in data:
                    if i[1]==atm and i[2]==fit and i[3]==rc and math.isnan(i[-1])==False:
                        err.append(i[-2])
                        err_abs.append(i[-1])
                        of.append(i)
                err_data.append([atm,fit,rc,round(mean(err),4),round(math.sqrt(variance(err)),4),round(mean(err_abs),4),round(math.sqrt(variance(err_abs)),4),round(median(err),4),round(median(err_abs),4)])
    x=[]
    y=[]
    yerr=[]
    c=[]
    yabs=[]
    yabs_rr=[]
    ymd=[]
    ymd_abs=[]
    for i in err_data:
        #x.append(f'{i[1]}-{i[2]}')
        x.append(i[2])
        y.append(i[3])
        yerr.append(i[4])
        yabs.append(i[5])
        yabs_rr.append(i[6])
        ymd.append(i[7])
        ymd_abs.append(i[8])
        c.append(i[1])
    width = 1200  # Logical width in pixels
    height = 800  # Logical height in pixels
    scale = 2

    fig = px.bar(x=x,y=y,error_y=yerr,color=c,barmode='group',
                 labels={'y':'Mean(Offset(RefDB)-Offset(LACS))','yerr':'Standard deviation of offset difference','x':'RC Model','color': 'Fit Model'})

    fig.write_html(f'{FIGS_PATH}/Mean_Offset_dif_barplot_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Mean_Offset_dif_barplot_{plot_atom}.jpeg',width=width, height=height,scale=scale)
    fig.write_image(f'{FIGS_PATH}/Mean_Offset_dif_barplot_{plot_atom}.pdf',width=width, height=height,scale=scale)


    fig = px.bar(x=x, y=yabs, error_y=yabs_rr, color=c, barmode='group',
                 labels={'y': 'Mean(|Offset(RefDB)-Offset(LACS)|)', 'yerr': 'Standard deviation of offset difference',
                         'x': 'RC Model', 'color': 'Fit Model'})
    fig.write_html(f'{FIGS_PATH}/Mean_Offset_abs_dif_barplot_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Mean_Offset_abs_dif_barplot_{plot_atom}.jpeg',width=width, height=height,scale=scale)
    fig.write_image(f'{FIGS_PATH}/Mean_Offset_abs_dif_barplot_{plot_atom}.pdf',width=width, height=height,scale=scale)


    fig = px.bar(x=x, y=ymd, error_y=yerr, color=c, barmode='group',
                 labels={'y': 'Median(Offset(RefDB)-Offset(LACS))', 'yerr': 'Standard deviation of offset difference',
                         'x': 'RC Model', 'color': 'Fit Model'})

    fig.write_html(f'{FIGS_PATH}/Median_Offset_dif_barplot_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Median_Offset_dif_barplot_{plot_atom}.jpeg', width=width, height=height, scale=scale)
    fig.write_image(f'{FIGS_PATH}/Median_Offset_dif_barplot_{plot_atom}.pdf', width=width, height=height, scale=scale)

    fig = px.bar(x=x, y=ymd_abs, error_y=yabs_rr, color=c, barmode='group',
                 labels={'y': 'Median(|Offset(RefDB)-Offset(LACS)|)', 'yerr': 'Standard deviation of offset difference',
                         'x': 'RC Model', 'color': 'Fit Model'})
    fig.write_html(f'{FIGS_PATH}/Median_Offset_abs_dif_barplot_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Median_Offset_abs_dif_barplot_{plot_atom}.jpeg', width=width, height=height, scale=scale)
    fig.write_image(f'{FIGS_PATH}/Median_Offset_abs_dif_barplot_{plot_atom}.pdf', width=width, height=height, scale=scale)


    offset_diff = [list(row) for row in zip(*of)]

    fig = px.histogram(x=offset_diff[-2],color=offset_diff[2],facet_col=offset_diff[3],facet_row=offset_diff[2],
                       labels={'x':'Offset(RefDB)-Offset(LACS)','facet_col':'RC Model','facet_row':'FIT Model'})
    fig.update_xaxes(title_font_size=10)
    fig.write_html(f'{FIGS_PATH}/Offset_dif_hist_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Offset_dif_hist_{plot_atom}.jpeg',width=width, height=height,scale=scale)
    fig.write_image(f'{FIGS_PATH}/Offset_dif_hist_{plot_atom}.pdf',width=width, height=height,scale=scale)


    fig = px.histogram(x=offset_diff[-1], color=offset_diff[2], facet_col=offset_diff[3], facet_row=offset_diff[2],
                        labels={'x': '|Offset(RefDB)-Offset(LACS)|', 'facet_col': 'RC Model', 'facet_row': 'FIT Model'})
    fig.update_xaxes(title_font_size=10)
    fig.write_html(f'{FIGS_PATH}/Offset_abs_dif_hist_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Offset_abs_dif_hist_{plot_atom}.jpeg',width=width, height=height,scale=scale)
    fig.write_image(f'{FIGS_PATH}/Offset_abs_dif_hist_{plot_atom}.pdf',width=width, height=height,scale=scale)

    fig = px.scatter(x=offset_diff[4], y=offset_diff[5],color=offset_diff[2], facet_col=offset_diff[3], facet_row=offset_diff[2],
                       labels={'x': 'Offset(RefDB)','y':'Offset(LACS)', 'facet_col': 'RC Model', 'facet_row': 'FIT Model'})
    #fig.update_xaxes(title_font_size=10)
    # for table s1
    # if plot_atom=='HA':
    #     xrange=0.4
    # elif plot_atom=='CA':
    #     xrange=3.5
    # elif plot_atom=='CB':
    #     xrange=3.0
    # elif plot_atom=='C':
    #     xrange=3.0
    # elif plot_atom=='N':
    #     xrange=8.5
    # else:
    #     print("ERROR")

    # for refdb all
    if plot_atom=='HA':
        xrange=0.6
    elif plot_atom=='CA':
        xrange=4.0
    elif plot_atom=='CB':
        xrange=4.0
    elif plot_atom=='C':
        xrange=6.0
    elif plot_atom=='N':
        xrange=10.0
    else:
        print("ERROR")


    fig.update_xaxes(range=[-xrange, xrange])
    fig.update_yaxes(range=[-xrange, xrange])
    fig.add_shape(
        type="line",
        x0=-xrange,  # Starting x-coordinate of the line
        y0=-xrange,  # Starting y-coordinate of the line
        x1=xrange,  # Ending x-coordinate of the line
        y1=xrange,  # Ending y-coordinate of the line
        line=dict(
            color="black",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        ),  row="all",
    col="all"
    )
    fig.write_html(f'{FIGS_PATH}/Offset_scatter_{plot_atom}.html')
    fig.write_image(f'{FIGS_PATH}/Offset_scatter_{plot_atom}.jpeg', width=width, height=height, scale=scale)
    fig.write_image(f'{FIGS_PATH}/Offset_scatter_{plot_atom}.pdf', width=width, height=height, scale=scale)

    # df = [list(row) for row in zip(*data)]
    # fig =px.scatter (x=df[4],y=df[5],color=df[1],facet_col=df[2],facet_row=df[3])
    # fig.update_xaxes(range=[-10,10])
    # fig.update_yaxes(range=[-10,10])
    # fig.show()
    # fig2 = px.histogram(x=df[-1],y=df[5],color=df[1],facet_col=df[2],facet_row=df[3],barmode='overlay')
    # fig2.show()


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
    #print (err,p)
    if data is not None:
        #print (data)
        for list_id in data:
            if fit in ['theilsen','ransac','quantile','tukey']:
                #for atom in data[list_id][f'offsets']:
                return data[list_id][f'offsets']
            else:
                #for atom in data[list_id][f'offsets_{fit}']:
                return data[list_id][f'offsets_{fit}']

def plot_cont_benchmark(csvfile):
    fit=[]
    rc=[]
    of_in=[]
    of_out_ca=[]
    of_out_cb = []
    of_out_c = []
    of_out_n = []
    dif_ca=[]
    dif_cb=[]
    dif_c=[]
    dif_n=[]
    abs_diff_ca=[]
    abs_diff_cb=[]
    abs_diff_c=[]
    abs_diff_n=[]
    with open(csvfile) as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            fit.append(row[0])
            rc.append(row[1])
            of_in.append(float(row[2]))
            of_out_ca.append(float(row[3]))
            of_out_cb.append(float(row[4]))
            of_out_c.append(float(row[5]))
            of_out_n.append(float(row[6]))
            dif_ca.append(float(row[2])-float(row[3]))
            dif_cb.append(float(row[2])-float(row[4]))
            dif_c.append(float(row[2])-float(row[5]))
            dif_n.append(float(row[2])-float(row[6]))
            abs_diff_ca.append(abs(float(row[2]) - float(row[3])))
            abs_diff_cb.append(abs(float(row[2]) - float(row[4])))
            abs_diff_c.append(abs(float(row[2]) - float(row[5])))
            abs_diff_n.append(abs(float(row[2]) - float(row[6])))

    fig = px.scatter(x=of_in,y=of_out_c,facet_col=fit,facet_row=rc,color=fit)
    xrange=4.0
    fig.update_xaxes(range=[-xrange, xrange])
    fig.update_yaxes(range=[-xrange, xrange])
    fig.add_shape(
        type="line",
        x0=-xrange,  # Starting x-coordinate of the line
        y0=-xrange,  # Starting y-coordinate of the line
        x1=xrange,  # Ending x-coordinate of the line
        y1=xrange,  # Ending y-coordinate of the line
        line=dict(
            color="black",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        ), row="all",
        col="all"
    )
    fig.show()
    fig = px.scatter(x=of_in,y=dif_c,facet_col=fit,facet_row=rc,color=fit)
    fig.show()



if __name__ == '__main__':
    #read_refdb_file('bmrpdbnew.txt')
    #plot_data('RefDB_benchmark.csv')
    #plot_data('RefDB_Lacs_benchmark.csv')
    plot_cont_benchmark('cbm.csv')