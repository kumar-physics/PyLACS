import numpy as np


class RandomCoil:


    def __init__(self):
        self.__three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                               'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                               'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                               'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
        self.__one_letter_code = dict([(value, key) for key, value in self.__three_letter_code.items()])

    def get(self, rc_name):
        if rc_name == 'wishart' or rc_name == 'wis':
            return self.wishart
        elif rc_name == 'wang' or rc_name == 'wan':
            return self.wang
        elif rc_name == 'lukhin' or rc_name == 'luk':
            return self.lukhin
        elif rc_name == 'schwarzinger' or rc_name == 'sch':
            return self.schwarzinger
        else:
            raise ValueError('Unknown random coil model')

    def get_value(self,res,atom,rc_name = None):
        res = res.upper()
        atom = atom.upper()
        if len(res)==3 and res not in self.__one_letter_code:
            raise ValueError(f'Unknown three letter code {res}')
        elif len(res)==1 and res not in self.__three_letter_code:
            raise ValueError(f'Unknown one letter code {res}')
        elif len(res) not in [1,3]:
            raise ValueError(f'Unknown residue {res}, please use one letter of three letter code')
        if res == 'GLY' and atom in ['HA2', 'HA3']:
            atom = 'HA'
        if atom not in self.atoms():
            raise ValueError(f'Unknown atom {atom} for the residue {res}')
        if rc_name is None:
            rc= self.get_average_all()
        elif type(rc_name) is list:
            rc = self.get_average(rc_name)
        else:
            rc = self.get(rc_name)
        if len(res) == 3:
            return rc[self.__one_letter_code[res]][self.atoms().index(atom)]
        else:
            return rc[res][self.atoms().index(atom)]


    def get_all(self):
        return self.wishart, self.wang, self.lukhin, self.schwarzinger

    def get_names(self, short_or_long = 'short'):
        if short_or_long == 'short':
            return ['wis','wan','luk','sch']
        elif short_or_long == 'long':
            return ['wishart', 'wang', 'lukhin', 'schwarzinger']
        else:
            raise ValueError('Unknown value for short_or_long')

    def atoms(self):
        return ['N','CO','CA','CB','H','HA']

    def get_average(self,name_list = None):
        if name_list is None:
            name_list = self.get_names()
        res = None
        sum = {}
        avg ={}
        # std = {}
        for name in name_list:
            rc_shifts = self.get(name)
            if res is None:
                residues = list(rc_shifts.keys())
            for res in residues:
                if res not in sum:
                    sum[res] = []
                sum[res].append([np.nan if item is None else item for item in rc_shifts[res]])
        for k in sum:
            avg[k] = [round(i,2) for i in np.mean(np.array(sum[k]), axis=0).tolist()]
            # std[k] = [round(i,2) for i in np.std(np.array(sum[k]), axis=0).tolist()]
        return avg

    def get_average_all(self):
        res = None
        sum = {}
        avg ={}
        # std = {}
        for name in self.get_names():
            rc_shifts = self.get(name)
            if res is None:
                residues = list(rc_shifts.keys())
            for res in residues:
                if res not in sum:
                    sum[res] = []
                sum[res].append([np.nan if item is None else item for item in rc_shifts[res]])
        for k in sum:
            avg[k] = [round(i,2) for i in np.mean(np.array(sum[k]), axis=0).tolist()]
            # std[k] = [round(i,2) for i in np.std(np.array(sum[k]), axis=0).tolist()]
        return avg

    def print(self):
        for name, rc_shifts in self.get_all().items():
            for res in rc_shifts:
                print (name, res, rc_shifts[res])


    wishart = {
                "I": [121.7, 176.4, 61.1, 38.8, 8.00, 4.17],
                "V": [120.5, 176.3, 62.2, 32.9, 8.03, 4.12],
                "D": [121.4, 176.3, 54.2, 41.1, 8.34, 4.64],
                "N": [119.0, 175.2, 53.1, 38.9, 8.40, 4.74],
                "F": [120.9, 175.8, 57.7, 39.6, 8.30, 4.62],
                "H": [118.2, 174.1, 55.0, 29.0, 8.42, 4.73],
                "W": [122.2, 176.1, 57.5, 29.60, 8.25, 4.66],
                "Y": [120.8, 175.9, 57.9, 38.8, 8.12, 4.55],
                "K": [121.6, 176.6, 56.2, 33.1, 8.29, 4.32],
                "L": [122.6, 177.6, 55.1, 42.4, 8.16, 4.34],
                "M": [120.7, 176.3, 55.4, 32.9, 8.28, 4.48],
                "Q": [120.6, 176.0, 55.7, 29.4, 8.32, 4.34],
                "R": [121.3, 176.3, 56.0, 30.9, 8.23, 4.34],
                "E": [121.7, 176.6, 56.6, 29.9, 8.42, 4.35],
                "T": [116.0, 174.7, 61.8, 69.8, 8.15, 4.35],
                "C": [118.9, 174.6, 58.2, 28.0, 8.43, 4.71],
                "S": [116.6, 174.6, 58.3, 63.8, 8.31, 4.47],
                "A": [125.0, 177.8, 52.5, 19.1, 8.24, 4.32],
                "G": [109.1, 174.9, 45.1, None, 8.33, 3.96],
                "P": [None, 177.3, 63.3, 32.1, None, 4.42],
                "B": [118.6, 174.6, 55.4, 41.1, 8.43, 4.71],
            }

    wang = {
                "I": [120.58, 175.52, 60.79, 38.43, 7.94, 4.18],
                "V": [119.91, 175.66, 62.00, 32.35, 7.98, 4.13],
                "D": [120.37, 176.00, 54.00, 40.78, 8.31, 4.62],
                "N": [118.50, 174.84, 53.00, 38.43, 8.35, 4.66],
                "F": [119.72, 175.46, 57.46, 39.41, 8.09, 4.59],
                "H": [118.92, 174.78, 55.74, 29.50, 8.18, 4.60],
                "W": [120.99, 175.87, 57.54, 29.60, 7.97, 4.60],
                "Y": [119.37, 175.29, 57.64, 38.78, 7.99, 4.56],
                "K": [121.10, 176.15, 56.29, 32.53, 8.17, 4.28],
                "L": [121.57, 176.70, 54.77, 42.14, 8.06, 4.36],
                "M": [120.14, 175.94, 55.43, 32.92, 8.22, 4.47],
                "Q": [119.82, 175.75, 55.89, 29.01, 8.20, 4.29],
                "R": [120.75, 176.01, 56.18, 30.36, 8.21, 4.26],
                "E": [120.62, 176.32, 56.66, 29.87, 8.36, 4.28],
                "T": [113.88, 174.78, 61.30, 68.92, 8.16, 4.44],
                "C": [118.10, 175.11, 58.24, 29.54, 8.10, 4.59],
                "S": [116.00, 174.41, 58.20, 63.75, 8.22, 4.45],
                "A": [123.82, 177.28, 52.46, 18.98, 8.09, 4.31],
                "G": [109.48, 174.01, 45.28, None, 8.37, 3.97],
                "P": [None, 176.62, 63.24, 31.81, None, 4.41],
                "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
            }

    lukhin = {
                "G": [110.03, 173.96, 45.41, None, 8.37, 3.97],
                "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
                "A": [125.10, 177.30, 52.42, 19.03, 8.09, 4.31],
                "S": [116.43, 174.41, 58.27, 64.14, 8.22, 4.45],
                "C": [119.19, 174.84, 58.01, 28.20, 8.43, 4.71],
                "M": [120.39, 175.45, 55.34, 33.00, 8.22, 4.47],
                "K": [121.82, 176.39, 56.59, 32.62, 8.21, 4.26],
                "V": [120.93, 175.79, 62.13, 32.65, 7.97, 4.60],
                "T": [114.17, 174.75, 61.62, 69.83, 8.16, 4.44],
                "I": [121.19, 175.69, 60.98, 38.87, 7.94, 4.18],
                "L": [122.22, 177.15, 54.82, 42.82, 8.06, 4.36],
                "D": [120.51, 176.45, 54.12, 40.83, 8.31, 4.62],
                "N": [119.17, 174.65, 53.22, 38.74, 8.35, 4.66],
                "E": [121.44, 176.27, 56.66, 30.13, 8.36, 4.28],
                "Q": [120.34, 175.54, 55.78, 29.34, 8.20, 4.29],
                "R": [122.42, 176.05, 56.25, 30.56, 8.21, 4.26],
                "H": [120.21, 174.54, 55.78, 29.78, 8.18, 4.60],
                "F": [119.84, 174.79, 57.91, 39.34, 8.09, 4.59],
                "Y": [120.41, 175.80, 57.77, 38.88, 7.99, 4.56],
                "W": [120.19, 175.85, 57.50, 29.09, 7.97, 4.60],
                "P": [136.73, 176.60, 63.27, 32.09, None, 4.41],
            }

    schwarzinger = {
                "A": [125.0, 178.5, 52.8, 19.3, 8.35, 4.35],
                "B": [118.7, 175.5, 55.6, 41.2, 8.54, 4.76],
                "C": [118.8, 175.3, 58.6, 28.3, 8.44, 4.59],
                "D": [119.1, 175.9, 53.0, 38.3, 8.56, 4.82],
                "E": [120.2, 176.8, 56.1, 29.9, 8.40, 4.42],
                "F": [120.7, 176.6, 58.1, 39.8, 8.31, 4.65],
                "G": [107.5, 174.9, 45.4, None, 8.41, 4.02],
                "H": [118.1, 175.1, 55.4, 29.1, 8.56, 4.79],
                "I": [120.4, 177.1, 61.6, 38.9, 8.17, 4.21],
                "K": [121.6, 177.4, 56.7, 33.2, 8.36, 4.36],
                "L": [122.4, 178.2, 55.5, 42.5, 8.28, 4.38],
                "M": [120.3, 177.1, 55.8, 32.9, 8.42, 4.52],
                "N": [119.0, 176.1, 53.3, 39.1, 8.51, 4.79],
                "P": [None, 177.8, 63.7, 32.2, None, 4.45],
                "Z": [None, None, 63.0, 34.8, None, 4.60],
                "Q": [120.5, 176.8, 56.2, 29.5, 8.44, 4.38],
                "R": [121.2, 177.1, 56.5, 30.9, 8.39, 4.38],
                "S": [115.5, 175.4, 58.7, 64.1, 8.43, 4.51],
                "T": [112.0, 175.6, 62.0, 70.0, 8.25, 4.43],
                "V": [119.3, 177.0, 62.6, 31.8, 8.16, 4.16],
                "W": [122.1, 177.1, 57.6, 29.8, 8.22, 4.70],
                "Y": [120.9, 176.7, 58.3, 38.9, 8.26, 4.58],
            }

if __name__ == "__main__":
    c = RandomCoil()
    print (c.get_value('G','H',['wis','wan']))

