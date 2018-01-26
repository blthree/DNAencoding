from utils import *
from functools import reduce
import logging

logging.basicConfig(level=logging.INFO)


dna_in = \
    ['ATAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTCG',
    'TATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAGC']
exp_dna = "TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGTACGTACGTACGTACGTACGACTAT"

exp_split_dna = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGA", "CATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTAC"]

exp_b3 = "20100202101010100021200012221110221201112000212210002212222212021100210122100110210111200021000000000000000000000000000010102"





def check_len_and_orientation(dna):
    valid_starts = ['A', 'T']
    valid_ends = ['G', 'C']
    for i in range(len(dna)):
        seq = dna[i]
        if len(seq) == 117:
            if not seq[0] in valid_starts and seq[-1] in valid_ends:
                seq_rc = rev_comp2(seq)
                if seq_rc[0] in valid_starts and seq_rc[-1] in valid_ends:
                    dna[i] = seq_rc[1:-1]
                else:
                    raise Exception("Invalid start/end chars")
            else:
                dna[i] = seq[1:-1]
        else:
            raise Exception("Sequence must be 117 nt")
    return dna


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


def dna_to_b3(str_to_decode, prev_char=None):
    dna_map = {'A': {'C': '0', 'G': '1', 'T': '2'},
               'C': {'G': '0', 'T': '1', 'A': '2'},
               'G': {'T': '0', 'A': '1', 'C': '2'},
               'T': {'A': '0', 'C': '1', 'G': '2'}}
    b3_out = ""
    for i in range(len(str_to_decode)):
        if not prev_char:
            prev_char = 'A'
        cur_char = dna_map[prev_char][str_to_decode[i]]
        b3_out += cur_char
        prev_char = str_to_decode[i]
    return b3_out

def b3_to_int(b3):
    out = 0
    for i in range(len(b3)):
        rev = b3[::-1]
        out += (3**i)*int(rev[i])
    return out


def confirm_parity(index):
    parity_trit = index[-1]
    file_id = index[:2]
    chunk_id = index[2:-1]
    if generate_parity_trit(chunk_id, file_id) == index:
        logging.debug("Parity Trit OK!")
    else:
        raise Exception("Parity Trit Error")
    out = (file_id, chunk_id)
    return out

def generate_parity_trit(i3, ID):
    i3 = "0"*(12-len(i3)) + i3
    temp = ID + i3
    even_trits = [int(temp[x]) for x in range(len(temp)) if x % 2 == 0]
    parity_trit = sum(even_trits) % 3
    index = ID + i3 + str(parity_trit)
    return index



def merge_overlapping(dna1, dna2):
    #print(len(dna1), len(dna2))
    #print("{0} \n{1}{2} ".format(dna1, "0"*(25+len(dna1)-len(dna2)) ,dna2))
    if dna1[-75:] == dna2[:75]:
        if len(dna1) % 25000 == 0:
            logging.info("Merged {0} sequences".format(((len(dna1)-100)/25)+1))
        return dna1 + dna2[-25:]
    else:
        raise Exception("Could not merge dna strings \n dna1: {0} \n dna2: {1} ".format(dna1[-75:], dna2[:75]))
def merge_overlapping2(dna1, dna2):
        if len(dna1) % 25000 == 0:
            logging.info("Merged {0} sequences".format(((len(dna1)-100)/25)+1))
        return dna1 + dna2[-25:]

def strip_len_info(b3):
    str_len = b3_to_int(b3[-20:])
    b3 = b3[:-20]
    zeros_to_remove = len(b3) - str_len
    assert b3[-zeros_to_remove:] == "0"*zeros_to_remove
    b3 = b3[:-zeros_to_remove]
    return b3


def b3_to_ord(b3):
    ords = bytearray()
    accum = ""
    for i in range(len(b3)):
        accum += b3[i]
        if accum in huffman_decode:
            ords.append(huffman_decode[accum])
            accum = ""
        else:
            pass
    return ords


dna_map = {'A': {'C': '0', 'G': '1', 'T': '2'},
           'C': {'G': '0', 'T': '1', 'A': '2'},
           'G': {'T': '0', 'A': '1', 'C': '2'},
           'T': {'A': '0', 'C': '1', 'G': '2'}}


with open("out.jpg", 'r') as f:
    metamorph = [line.strip("\n") for line in f.readlines()]


####### BEGIN decode metamorph

verified_dna = check_len_and_orientation(metamorph)
indexed_dna = dict()
indexed_b3 = dict()
for seq in verified_dna:
    indexed_dna[seq[-15:]] = seq[:-15]
    indexed_b3[dna_to_b3(seq[-15:], prev_char=seq[-16])] = seq[:-15]

parsed_index = {}
for x in indexed_b3:
    parsed_index[confirm_parity(x)] = indexed_b3[x]
files = dict()
for key in parsed_index:
    file_id, chunk_id = key
    if file_id not in files:
        files[file_id] = {b3_to_int(chunk_id): parsed_index[key]}
    else:
        files[file_id][b3_to_int(chunk_id)] = parsed_index[key]
f1 = files["12"]
for key in f1:
    if key % 2 != 0:
        f1[key] = rev_comp2(f1[key])

def split_up(dna_list):
    groups = list()
    ordered = [dna_list[i] for i in range(len(dna_list))]
    while len(ordered) > 0:
        if len(ordered) < 10000:
            groups.append(ordered)
            ordered = list()
        else:
            groups.append(ordered[:10000])
            ordered = ordered[10000:]
    return groups

logging.info("Found {0} sequences, encoding approximately {1} characters".format(len(f1), len(f1)*18))
# TODO: make me faster, need to dump to disk or something
# if we just split and map reduce, then we lose the previous char!
# could also convert then merge, by passing a tuple of (current seq, prev_seq[:-1])
split_up(f1)
merged = reduce(merge_overlapping2, [f1[k] for k in range(len(f1))])

logging.info("Merge complete")
logging.info("Begin dna to base3 conversion")
b3_message = dna_to_b3(merged)

ord_data = b3_to_ord(strip_len_info(b3_message))

#decoded2 = map(chr, ord_data)

with open("decoded.jpg", "wb") as f:
    f.write(ord_data)

