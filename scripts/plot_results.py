import plotly.express as px
import csv
from pandas import DataFrame

def plot_csv(csvfile):
    with open(csvfile,'r') as csvFile:
        csvf = csv.reader(csvFile)
        data = list(csvf)
    clean_data = [x for x in data if len(x)==31]
    for x in data:
        print (len(x))
        if len(x)!=31:
            print (x)
    result = [list(x) for x in zip(*clean_data)]
    print (result)
    df = DataFrame({#["tukey", "theilsen", "ransac", "quantile", "bayes"]
        'BMRB_ID':result[0],
        'RefDB(CA)':[float(i) if i!='' else None for i in result[1]],
        'RefDB(C)':[float(i) if i!='' else None for i in result[2]],
        'RefDB(HA)':[float(i) if i!='' else None for i in result[3]],
        'LACS(CA)': [float(i) if i != '' else None for i in result[4]],
        'LACS(C)': [float(i) if i != '' else None for i in result[5]],
        'LACS(HA)': [float(i) if i != '' else None for i in result[6]],
        'Tukey(CA)': [float(i) if i!='' else None for i in result[7]],
        'Tukey(C)': [float(i) if i!='' else None for i in result[8]],
        'Tukey(HA)': [float(i) if i!='' else None for i in result[9]],
        'Theilsen(CA)': [float(i) if i != '' else None for i in result[10]],
        'Theilsen(C)': [float(i) if i != '' else None for i in result[11]],
        'Theilsen(HA)': [float(i) if i != '' else None for i in result[12]],
        'Ransac(CA)': [float(i) if i != '' else None for i in result[13]],
        'Ransac(C)': [float(i) if i != '' else None for i in result[14]],
        'Ransac(HA)': [float(i) if i != '' else None for i in result[15]],
        'Quantile(CA)': [float(i) if i != '' else None for i in result[16]],
        'Quantile(C)': [float(i) if i != '' else None for i in result[17]],
        'Quantile(HA)': [float(i) if i != '' else None for i in result[18]],
        'Bayes(CA)': [float(i) if i != '' else None for i in result[19]],
        'Bayes(C)': [float(i) if i != '' else None for i in result[20]],
        'Bayes(HA)': [float(i) if i != '' else None for i in result[21]],
        'Bayes(CA_SD)': [float(i) if i != '' else None for i in result[24]],
        'Bayes(C_SD)': [float(i) if i != '' else None for i in result[27]],
        'Bayes(HA_SD)': [float(i) if i != '' else None for i in result[30]],



    })
    fig=px.scatter(df,x='RefDB(CA)',y=['LACS(CA)','Tukey(CA)','Theilsen(CA)','Ransac(CA)','Quantile(CA)','Bayes(CA)'],hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-2,  # Starting x-coordinate of the line
        y0=-2,  # Starting y-coordinate of the line
        x1=3,  # Ending x-coordinate of the line
        y1=3,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_CA_Scatter.html')
    fig = px.scatter(df, x='RefDB(HA)',
                     y=['LACS(HA)', 'Tukey(HA)', 'Theilsen(HA)', 'Ransac(HA)', 'Quantile(HA)', 'Bayes(HA)'],
                     hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-0.5,  # Starting x-coordinate of the line
        y0=-0.5,  # Starting y-coordinate of the line
        x1=0.4,  # Ending x-coordinate of the line
        y1=0.4,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_HA_Scatter.html')
    fig = px.scatter(df, x='RefDB(C)',
                     y=['LACS(C)', 'Tukey(C)', 'Theilsen(C)', 'Ransac(C)', 'Quantile(C)', 'Bayes(C)'],
                     hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-3,  # Starting x-coordinate of the line
        y0=-3,  # Starting y-coordinate of the line
        x1=4,  # Ending x-coordinate of the line
        y1=4,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_C_Scatter.html')
    fig = px.scatter(df, x='RefDB(CA)', y = 'Bayes(CA)', error_y='Bayes(CA_SD)',hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-2,  # Starting x-coordinate of the line
        y0=-2,  # Starting y-coordinate of the line
        x1=3,  # Ending x-coordinate of the line
        y1=3,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_Bayes_CA_Scatter.html')
    fig = px.scatter(df, x='RefDB(C)', y='Bayes(C)', error_y='Bayes(C_SD)', hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-3,  # Starting x-coordinate of the line
        y0=-3,  # Starting y-coordinate of the line
        x1=4,  # Ending x-coordinate of the line
        y1=4,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_Bayes_C_Scatter.html')
    fig = px.scatter(df, x='RefDB(HA)', y='Bayes(HA)', error_y='Bayes(HA_SD)', hover_name='BMRB_ID')
    fig.add_shape(
        type="line",
        x0=-0.5,  # Starting x-coordinate of the line
        y0=-0.5,  # Starting y-coordinate of the line
        x1=0.4,  # Ending x-coordinate of the line
        y1=0.4,  # Ending y-coordinate of the line
        line=dict(
            color="red",
            width=2,
            dash="dash"  # Optional: 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'
        )
    )
    fig.write_html('RefDBvsLACS_Bayes_HA_Scatter.html')

if __name__ == "__main__":
    plot_csv('results2.csv')