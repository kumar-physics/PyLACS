import math
import csv
from typing import Mapping
import plotly.express as px
import pandas as pd
from pandas.core.interchange.dataframe_protocol import DataFrame

from random_coil import RandomCoil
from sklearn.neighbors import KernelDensity
import numpy as np
from scipy.stats import gaussian_kde

_THREE_TO_ONE: Mapping[str, str] = {
        'ILE': 'I', 'GLN': 'Q', 'GLY': 'G', 'GLU': 'E', 'CYS': 'C',
        'ASP': 'D', 'SER': 'S', 'LYS': 'K', 'PRO': 'P', 'ASN': 'N',
        'VAL': 'V', 'THR': 'T', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
        'ALA': 'A', 'MET': 'M', 'LEU': 'L', 'ARG': 'R', 'TYR': 'Y',
    }
_ONE_TO_THREE: Mapping[str, str] = {v: k for k, v in _THREE_TO_ONE.items()}

def read_csv(fname):
    rc=RandomCoil()
    res={}
    atm={}
    with open(fname,'r') as csvFile:
        csv_data = csv.reader(csvFile)
        for row in csv_data:
            if row[4] in _THREE_TO_ONE:
                if f'{row[4]}-{row[6]}' not in res:
                    res[f'{row[4]}-{row[6]}']=[]
                if row[6] not in atm:
                    atm[row[6]]=[]
                # res[f'{row[4]}-{row[6]}'].append(rc.get_value(row[4],row[6])-float(row[-1]))
                # atm[row[6]].append(rc.get_value(row[4],row[6])-float(row[-1]))
                res[f'{row[4]}-{row[6]}'].append(float(row[-1]))
                atm[row[6]].append(float(row[-1]))
    return res, atm


def cal_entropy(bf,af):
    bf_res, bf_atm = read_csv(bf)
    af_res, af_atm = read_csv(af)
    r=[]
    a=[]
    sb=[]
    sa=[]
    data={}

    data1={'Res-Atom':[],'Entropy':[]}
    data2={'Res-Atom':[],'Entropy':[]}
    data3={'Res-Atom':[],'Entropy':[]}
    for k in bf_res:
        # print (k,entropy(bf_res[k]),entropy(af_res[k]))
        # print(k, differential_entropy_kde(bf_res[k]), differential_entropy_kde(af_res[k]))

        data1['Res-Atom'].append(k)
        data2['Res-Atom'].append(k)
        data3['Res-Atom'].append(k)
        eb=differential_entropy_kde(bf_res[k])
        ea=differential_entropy_kde(af_res[k])
        data1['Entropy'].append(eb)
        data2['Entropy'].append(ea)
        data3['Entropy`'].append(eb-ea)
        d1 = {'Chemical Shift [ppm]': bf_res[k]}
        d2 = {'Chemical Shift [ppm]': af_res[k]}
        f1 = pd.DataFrame(d1)
        f1['Dataset'] = f'Before (Entropy ={round(eb,3)})'
        f2 = pd.DataFrame(d2)
        f2['Dataset'] = f'After (Entropy ={round(ea,3)})'
        df = pd.concat([f1, f2])
        fig1 = px.histogram(df, x='Chemical Shift [ppm]', color='Dataset',title=k,barmode='overlay')
        fig1.update_xaxes(autorange='reversed')
        fig1.write_html(f'/Users/kumaranbaskaran/entropy/entropy_{k}.html')
    df1 = pd.DataFrame(data1)
    df1['Dataset'] = 'Before'
    df2 = pd.DataFrame(data2)
    df2['Dataset'] = 'After'
    df3 = pd.DataFrame(data3)
    df3['Dataset'] = 'Before-After'
    all_df = pd.concat([df1,df2,df3])
    df_sorted = all_df.sort_values(by='Res-Atom', ascending=True)
    fig = px.bar(df_sorted,x='Res-Atom',y='Entropy',color = 'Dataset',barmode='group')
    fig.write_html('/Users/kumaranbaskaran/entropy/Entropy_comparison1.html')
    data1 = {'Atom': [], 'Entropy': []}
    data2 = {'Atom': [], 'Entropy': []}
    for k in bf_atm:
        # print (k,entropy(bf_atm[k]),entropy(af_atm[k]))
        # print(k, differential_entropy_kde(bf_atm[k]), differential_entropy_kde(af_atm[k]))
        data1['Atom'].append(k)
        data2['Atom'].append(k)
        data1['Entropy'].append(differential_entropy_kde(bf_atm[k]))
        data2['Entropy'].append(differential_entropy_kde(af_atm[k]))
    df1 = pd.DataFrame(data1)
    df1['Dataset'] = 'Before'
    df2 = pd.DataFrame(data2)
    df2['Dataset'] = 'After'
    all_df = pd.concat([df1, df2])

    fig = px.bar(all_df, x='Atom', y='Entropy', color='Dataset', barmode='group')
    fig.write_html('/Users/kumaranbaskaran/entropy/Entropy_comparison2.html')



def entropy(dist):
    """
    Calculate the Shannon entropy of a distribution.

    Parameters
    ----------
    dist : list of float
        Distribution values (can be counts/frequencies or probabilities).

    Returns
    -------
    float
        Entropy in bits (log base 2).
    """
    total = sum(dist)
    if total == 0:
        return 0.0

    # Normalize to probabilities
    probs = [x / total for x in dist if x > 0]

    # Compute entropy
    return round(-sum(p * math.log2(p) for p in probs),5)






def _clean(samples):
    x = np.asarray(samples, dtype=float)
    x = x[np.isfinite(x)]            # drop NaN/±inf
    return x

def discrete_entropy(samples, bins=30, base=2):
    """
    Shannon entropy (discrete) estimated via histogram.
    Robust to zeros and degenerate inputs.
    """
    x = _clean(samples)
    if x.size == 0:
        return np.nan  # no data

    counts, _ = np.histogram(x, bins=bins, density=False)
    total = counts.sum()
    if total == 0:
        return np.nan  # all bins empty (e.g., all NaN originally)

    p = counts / total
    mask = p > 0
    if not np.any(mask):
        return 0.0

    H = -np.sum(p[mask] * np.log(p[mask])) / np.log(base)
    return H

def differential_entropy_kde(samples, base=2, grid=2048, eps=1e-300, jitter=0.0):
    """
    Differential entropy (continuous) via KDE and Riemann sum:
        H = -∫ p(x) log p(x) dx
    Handles constant-data cases and renormalizes the KDE on the grid.
    """
    x = _clean(samples)
    if x.size == 0:
        return np.nan

    # If all samples identical, gaussian_kde will fail (singular covariance).
    if np.allclose(x, x[0]):
        if jitter > 0:
            x = x + np.random.normal(scale=jitter, size=x.size)
        else:
            # For a delta distribution, true differential entropy → -inf.
            # Return a large negative number instead of -inf to avoid explosions.
            return -1e9

    kde = gaussian_kde(x)

    # Choose a grid wide enough (±3σ) to capture tails
    mu = np.mean(x)
    sd = np.std(x, ddof=1)
    if sd == 0:  # safety (should be caught above)
        sd = 1.0

    lo = mu - 3*sd
    hi = mu + 3*sd
    xs = np.linspace(lo, hi, grid)
    px = kde(xs)

    # Clip to avoid exact zeros, then renormalize so ∑ px * dx ≈ 1
    px = np.clip(px, eps, None)
    dx = xs[1] - xs[0]
    Z = np.sum(px) * dx
    if Z <= 0 or not np.isfinite(Z):
        return np.nan
    px /= Z

    H = -np.sum(px * np.log(px)) * dx / np.log(base)
    return H



if __name__ == "__main__":
    cal_entropy('/Users/kumaranbaskaran/entropy/cs_before.csv','/Users/kumaranbaskaran/entropy/cs_after.csv')