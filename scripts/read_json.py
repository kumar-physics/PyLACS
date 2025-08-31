import csv
import json

def read_csv_file(csvfile):
    fo=open('results2.csv','w')
    with open(csvfile,'r') as csvFile:
        csvF = csv.reader(csvFile)
        for row in csvF:
            bmrb_id = row[0]
            print (bmrb_id)
            l=f'{bmrb_id},{row[4]},{row[6]},{row[8]},{row[5]},{row[7]},{row[9]}'
            if bmrb_id not in  ['4077','5168']:
                for method in ["tukey", "theilsen", "ransac", "quantile", "bayes"]:
                    with open(f'../scratch/scratch/{bmrb_id}/{bmrb_id}_{method}.json','r') as f:
                        data = json.load(f)
                        of = data['1']['offsets']
                        try:
                            ca = -of['ca']
                        except KeyError:
                            cs = ''
                        try:
                            c = -of['c']
                        except KeyError:
                            c=''
                        try:
                            h= -of['ha']
                        except KeyError:
                            h=''
                        l = f'{l},{ca},{c},{h}'
                        if method == 'bayes':
                            ofb = data['1']['offsets_bayes']
                            try:
                                ca_min = ofb['ca']['ci95'][0]
                                ca_max = ofb['ca']['ci95'][1]
                                ca_sd = ofb['ca']['sd']
                            except KeyError:
                                ca_min=''
                                ca_max = ''
                                ca_sd = ''
                            try:
                                c_min = ofb['c']['ci95'][0]
                                c_max = ofb['c']['ci95'][1]
                                c_sd = ofb['c']['sd']
                            except KeyError:
                                c_min=''
                                c_max = ''
                                c_sd = ''
                            try:
                                ha_min = ofb['ha']['ci95'][0]
                                ha_max = ofb['ha']['ci95'][1]
                                ha_sd = ofb['ha']['sd']
                            except KeyError:
                                ha_min=''
                                ha_max = ''
                                ha_sd = ''
                            l = f'{l},{ca_min},{ca_max},{ca_sd},{c_min},{c_max},{c_sd},{ha_min},{ha_max},{ha_sd}'


                fo.write(f'{l}\n')



if __name__ == "__main__":
    read_csv_file('Table_S1_clean.csv')