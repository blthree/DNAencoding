
class asciiStr(object):
    huffman = {0: '22201', 85: '22200', 170: '22122', 127: '22121', 253: '22120', 52: '22112', 138: '22111',
               41: "22110",
               86: '22102', 42: '22101', 100: '22100', 44: '22022', 250: '22020', 132: '22021', 161: '22012',
               98: '22010',
               8: '22002', 34: '22011', 10: '22001', 149: '22000', 87: '21222', 21: '21221', 74: '21220', 36: '21212',
               69: '21210', 177: '21202', 20: '21211', 213: '21200', 163: '21201', 229: '21121', 255: '21122',
               197: '21120',
               133: '21112', 252: '21110', 26: '21111', 173: '21101', 151: '21102', 82: '21100', 75: '21022',
               37: '21021',
               166: '21011', 191: '21020', 88: '21012', 63: '21010', 68: '21001', 150: '21002', 76: '21000', 4: '20222',
               154: '20221', 234: '20212', 22: '20220', 162: '20211', 105: '20210', 102: '20202', 171: '20201',
               104: '20200', 169: '20122', 196: '20121', 208: '20120', 84: '20112', 130: '20111', 146: '20102',
               72: '20110',
               16: '20101', 66: '20100', 24: '20022', 106: '20012', 223: '20020', 58: '20021', 137: '20011',
               73: '20010',
               101: '20001', 168: '20002', 181: '12221', 175: '12222', 251: '20000', 40: '12220', 140: '12212',
               17: '12211',
               83: '12210', 254: '12202', 240: '12201', 214: '12200', 53: '12122', 202: '12112', 25: '12121',
               18: '12120',
               247: '12111', 174: '12110', 112: '12102', 89: '12101', 210: '12100', 217: '12012', 248: '12020',
               194: '12021', 182: '12022', 80: '12011', 79: '12002', 195: '12010', 12: '12001', 209: '12000',
               165: '11222',
               245: '11221', 2: '11220', 81: '11212', 38: '11211', 141: '11202', 211: '11210', 239: '11200',
               95: '11201',
               43: '11122', 224: '11121', 203: '11112', 145: '11120', 147: '11110', 19: '11111', 50: '11101',
               136: '11102',
               107: '11100', 134: '11022', 109: '11021', 153: '11020', 148: '11002', 205: '11010', 212: '11011',
               54: '11012', 241: '11000', 156: '11001', 115: '10222', 116: '10221', 78: '10220', 67: '10211',
               70: '10212',
               178: '10210', 159: '10202', 142: '10201', 92: '10200', 48: '10122', 90: '10120', 218: '10121',
               126: '10112',
               39: '10111', 219: '10102', 167: '10110', 114: '10101', 172: '10022', 14: '10100', 120: '10020',
               139: '10021',
               160: '10012', 33: '10011', 179: '10010', 117: '10002', 225: '10001', 129: '10000', 183: '02222',
               230: '02220', 35: '02221', 93: '02210', 6: '02211', 32: '02212', 56: '02201', 158: '02202', 185: '02121',
               47: '02122', 143: '02200', 123: '02111', 204: '02120', 242: '02112', 111: '02110', 103: '02102',
               108: '02101', 9: '02100', 65: '02022', 249: '02020', 13: '02021', 180: '02012', 226: '02001',
               144: '02002',
               15: '02010', 57: '02011', 128: '02000', 135: '01220', 243: '01221', 190: '01222', 207: '01212',
               77: '01211',
               45: '01210', 91: '01202', 192: '01201', 186: '01122', 216: '01200', 97: '01112', 118: '01120',
               246: '01121',
               215: '01111', 51: '01102', 206: '01110', 184: '01100', 227: '01101', 233: '01022', 237: '01021',
               188: '01020', 113: '01012', 49: '01011', 201: '01010', 155: '01002', 222: '01000', 231: '01001',
               5: '00222',
               27: '00221', 131: '00212', 164: '00220', 3: '00211', 46: '00210', 119: '00201', 28: '00202',
               176: '00200',
               23: '00122', 64: '00121', 157: '00120', 187: '00112', 244: '00110', 238: '00111', 96: '00102',
               235: '00101',
               60: '00022', 1: '00100', 110: '00021', 200: '00011', 221: '00020', 99: '00012', 31: '00010',
               198: '00002',
               193: '00001', 125: '00000', 124: '222222', 152: '222221', 122: '222220', 71: '222212', 94: '222211',
               220: '222210', 29: '222202', 199: '222201', 61: '222200', 11: '222122', 228: '222121', 62: '222120',
               55: '222112', 121: '222111', 7: '222110', 30: '222102', 232: '222101', 189: '222100', 59: '222021',
               236: '222022'}
    dna_map = {'A': {'0': 'C', '1': 'G', '2': 'T'},
               'C': {'0': 'G', '1': 'T', '2': 'A'},
               'G': {'0': 'T', '1': 'A', '2': 'C'},
               'T': {'0': 'A', '1': 'C', '2': 'G'}}
    def __init__(self, in_str):
        if not isinstance(in_str, bytes):
            self.asciiChars = [ord(l) for l in in_str]
        else:
            self.asciiChars = [int(x) for x in in_str]

        self.raw_trits = ""
        self.file_idx = ""
        self.dna = ""
        return

    def _a2b3(self):
        for ascii_code in self.asciiChars:
            yield self.huffman[ascii_code]

    def _add_zeros(self, mod_num=25):
        b3_str_len = len(self.raw_trits)
        b3_len_enc = self._int_to_b3(b3_str_len)
        b3_len_enc = "0" * (20 - len(b3_len_enc)) + b3_len_enc
        zeros_to_add = ((((b3_str_len + len(b3_len_enc)) // mod_num) + 1) * mod_num - (b3_str_len + len(b3_len_enc)))
        self.file_idx = "0" * zeros_to_add + b3_len_enc
        self.trits = self.raw_trits + self.file_idx
        return None

    def b3_to_dna(self, prev_char=None):
        dna_out = ""
        for i in range(len(self.trits)):
            if not prev_char:
                prev_char = 'A'
            cur_char = self.dna_map[prev_char][self.trits[i]]
            dna_out += cur_char
            prev_char = cur_char
        self.dna = dna_out


    @staticmethod
    def from_file(filename):
        with open(filename, 'rb') as f:
            return asciiStr(f.read())

    def _int_to_b3(self, i):
        out = ""
        while i > 0:
            rem = i % 3
            out = str(rem) + out
            i = i // 3
        return out

    def to_trits(self):
        self.raw_trits = "".join(self._a2b3())
        self._add_zeros()
        self.b3_to_dna()
        # next step is split and index
        pass

class indexed_dna(object):
    def __init__(self, seq, idx):
        if len(seq) == 125:
            # split off idx
            self.idx_dna = idx

        pass

class trits(object):
    def __init__(self, trits):
        self.trits = trits


class dna(object):
    complements = str.maketrans({'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'})
    dna_map = {'A': {'C': '0', 'G': '1', 'T': '2'},
               'C': {'G': '0', 'T': '1', 'A': '2'},
               'G': {'T': '0', 'A': '1', 'C': '2'},
               'T': {'A': '0', 'C': '1', 'G': '2'}}

    def __init__(self, seq):
        self.seq = seq

    def rev_comp(self):
        self.seq = self.seq[::-1].translate(self.complements)
        return None

    def to_base3(self, prev_char=None):
        b3_out = ""
        for i in range(len(self.seq)):
            if not prev_char:
                prev_char = 'A'
            b3_out += self.dna_map[prev_char][self.seq[i]]
            prev_char = self.seq[i]
        return b3_out