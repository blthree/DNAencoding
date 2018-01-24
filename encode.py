import binascii
from random import randint

s0 = "Birney and Goldman"
exp_dna_out = "TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGC" \
          "TGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGTACGTACGTACGTACGTACGACTAT"
exp_split_dna = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGA", "CATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTAC"]

exp_115 = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTC",  "ATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAG"]


huffman = {0: '22201', 85: '22200', 170: '22122', 127: '22121', 253: '22120', 52: '22112', 138: '22111', 41: "22110",
           86: '22102', 42: '22101', 100: '22100', 44: '22022', 250: '22020', 132: '22021', 161: '22012', 98: '22010',
           8: '22002', 34: '22011', 10: '22001', 149: '22000', 87: '21222', 21: '21221', 74: '21220', 36: '21212',
           69: '21210', 177: '21202', 20: '21211', 213: '21200', 163: '21201', 229: '21121', 255: '21122', 197: '21120',
           133: '21112', 252: '21110', 26: '21111', 173: '21101', 151: '21102', 82: '21100', 75: '21022', 37: '21021',
           166: '21011', 191: '21020', 88: '21012', 63: '21010', 68: '21001', 150: '21002', 76: '21000', 4: '20222',
           154: '20221', 234: '20212', 22: '20220', 162: '20211', 105: '20210', 102: '20202', 171: '20201',
           104: '20200', 169: '20122', 196: '20121', 208: '20120', 84: '20112', 130: '20111', 146: '20102', 72: '20110',
           16: '20101', 66: '20100', 24: '20022', 106: '20012', 223: '20020', 58: '20021', 137: '20011', 73: '20010',
           101: '20001', 168: '20002', 181: '12221', 175: '12222', 251: '20000', 40: '12220', 140: '12212', 17: '12211',
           83: '12210', 254: '12202', 240: '12201', 214: '12200', 53: '12122', 202: '12112', 25: '12121', 18: '12120',
           247: '12111', 174: '12110', 112: '12102', 89: '12101', 210: '12100', 217: '12012', 248: '12020',
           194: '12021', 182: '12022', 80: '12011', 79: '12002', 195: '12010', 12: '12001', 209: '12000', 165: '11222',
           245: '11221', 2: '11220', 81: '11212', 38: '11211', 141: '11202', 211: '11210', 239: '11200', 95: '11201',
           43: '11122', 224: '11121', 203: '11112', 145: '11120', 147: '11110', 19: '11111', 50: '11101', 136: '11102',
           107: '11100', 134: '11022', 109: '11021', 153: '11020', 148: '11002', 205: '11010', 212: '11011',
           54: '11012', 241: '11000', 156: '11001', 115: '10222', 116: '10221', 78: '10220', 67: '10211', 70: '10212',
           178: '10210', 159: '10202', 142: '10201', 92: '10200', 48: '10122', 90: '10120', 218: '10121', 126: '10112',
           39: '10111', 219: '10102', 167: '10110', 114: '10101', 172: '10022', 14: '10100', 120: '10020', 139: '10021',
           160: '10012', 33: '10011', 179: '10010', 117: '10002', 225: '10001', 129: '10000', 183: '02222',
           230: '02220', 35: '02221', 93: '02210', 6: '02211', 32: '02212', 56: '02201', 158: '02202', 185: '02121',
           47: '02122', 143: '02200', 123: '02111', 204: '02120', 242: '02112', 111: '02110', 103: '02102',
           108: '02101', 9: '02100', 65: '02022', 249: '02020', 2021: '2136', 180: '02012', 226: '02001', 144: '02002',
           15: '02010', 57: '02011', 128: '02000', 135: '01220', 243: '01221', 190: '01222', 207: '01212', 77: '01211',
           45: '01210', 91: '01202', 192: '01201', 186: '01122', 216: '01200', 97: '01112', 118: '01120', 246: '01121',
           215: '01111', 51: '01102', 206: '01110', 184: '01100', 227: '01101', 233: '01022', 237: '01021',
           188: '01020', 113: '01012', 49: '01011', 201: '01010', 155: '01002', 222: '01000', 231: '01001', 5: '00222',
           27: '00221', 131: '00212', 164: '00220', 3: '00211', 46: '00210', 119: '00201', 28: '00202', 176: '00200',
           23: '00122', 64: '00121', 157: '00120', 187: '00112', 244: '00110', 238: '00111', 96: '00102', 235: '00101',
           60: '00022', 1: '00100', 110: '00021', 200: '00011', 221: '00020', 99: '00012', 31: '00010', 198: '00002',
           193: '00001', 125: '00000', 124: '222222', 152: '222221', 122: '222220', 71: '222212', 94: '222211',
           220: '222210', 29: '222202', 199: '222201', 61: '222200', 11: '222122', 228: '222121', 62: '222120',
           55: '222112', 121: '222111', 7: '222110', 30: '222102', 232: '222101', 189: '222100', 59: '222021',
           236: '222022'}


def a2b3(s):
    for ascii_code in s:
        yield huffman[ascii_code]

def int_to_b3(i):
    out = ""
    while i > 0:
        rem = i % 3
        out = str(rem) + out
        i = i // 3
    return out

def build_initial_dna_string(b3_str):
    # calculate length for indexing
    str_to_encode = add_zeros(b3_str)

    return b3_to_dna(str_to_encode)

def b3_to_dna(str_to_encode, prev_char=None):

    dna_map = {'A': {'0': 'C', '1': 'G', '2': 'T'},
               'C': {'0': 'G', '1': 'T', '2': 'A'},
               'G': {'0': 'T', '1': 'A', '2': 'C'},
               'T': {'0': 'A', '1': 'C', '2': 'G'}}
    dna_out = ""
    for i in range(len(str_to_encode)):
        if not prev_char:
            prev_char = 'A'
        cur_char = dna_map[prev_char][str_to_encode[i]]
        dna_out += cur_char
        prev_char = cur_char

    return dna_out

def add_zeros(b3_str, mod_num=25):
    b3_str_len = len(b3_str)
    b3_len_enc = int_to_b3(b3_str_len)
    b3_len_enc = "0" * (20 - len(b3_len_enc)) + b3_len_enc
    zeros_to_add = ((((b3_str_len + len(b3_len_enc)) // mod_num) + 1) * mod_num - (b3_str_len + len(b3_len_enc)))
    str_to_encode = b3_str + "0" * zeros_to_add + b3_len_enc
    assert (len(str_to_encode) % mod_num) == 0
    return str_to_encode

def chunk_str(s, max_len=100, o=25):
    """
    chunks a string into a list of strings of at most max_len
        with an overlap of o letters
    :param s: string to split
    :param max_len: length of each string
    :param o: overlap between chunks
    :return: list of strings
    """
    start = 0
    end = max_len
    out = list()
    while start <= len(s)-max_len:
        out.append(s[start:start+end])
        start += o
    return out

def rev_comp(dna):
    complements = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    dna_out = ""
    for l in dna[::-1]:
        dna_out += complements[l]
    return dna_out

def split_and_index(dna):
    dna = chunk_str(dna)
    for i in range(len(dna)):
        if i % 2 == 1:
            dna[i] = rev_comp(dna[i])
    return dna

def generate_parity_trit(dna, i, ID):
    i3 = int_to_b3(i)
    i3 = "0"*(12-len(i3)) + i3
    temp = ID + i3
    even_trits = [int(temp[x]) for x in range(len(temp)) if x % 2 == 0]
    parity_trit = sum(even_trits) % 3
    index = ID + i3 + str(parity_trit)
    return index

def add_final_trits(dna):
    if dna.startswith("A"):
        dna = "T" + dna
    elif dna.startswith("T"):
        dna = "A" + dna
    else:
        dna = ["A", "T"][randint(0, 10) % 2]
    if dna.endswith("C"):
        dna += "G"
    elif dna.endswith("G"):
        dna += "C"
    else:
        dna += ["C", "G"][randint(0, 10) % 2]
    return dna

s0_ord = [ord(l) for l in s0]
s1 = "".join(a2b3(s0_ord))
s2 = build_initial_dna_string(s1)
assert s2 == exp_dna_out
print(s2)
print(exp_dna_out)
s5 = split_and_index(s2)
assert rev_comp("AATTGCA") == "TGCAATT"
ID = "12"
for i in range(len(s5)):
    s5[i] += b3_to_dna(generate_parity_trit(s5[i], i, ID), prev_char=s5[i][-1])
    assert len(s5[i]) == 115
print(s5)
print(exp_115)
assert s5 == exp_115
for i in range(len(s5)):
    s5[i] = add_final_trits(s5[i])
    assert len(s5[i]) == 117
# this is the final 117 bp dna string
# will be 125 bp with illumina adapters
print(s5)
