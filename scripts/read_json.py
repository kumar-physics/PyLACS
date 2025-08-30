import csv
import json

def read_csv_file(csvfile):
    fo=open('results.csv','w')
    with open(csvfile,'r') as csvFile:
        csvF = csv.reader(csvFile)
        for row in csvF:
            bmrb_id = row[0]
            print (bmrb_id)
            l=f'{bmrb_id},{row[4]},{row[6]},{row[8]},{row[5]},{row[7]},{row[9]}'
            if bmrb_id not in  ['4077','5168']:
                for method in ["tukey", "theilsen", "ransac", "quantile", "bayes"]:
                    with open(f'../scratch/{bmrb_id}/{bmrb_id}_{method}.json','r') as f:
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
            fo.write(f'{l}\n')



if __name__ == "__main__":
    read_csv_file('Table_S1_clean.csv')