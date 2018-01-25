import binascii
from random import randint
from utils import huffman

s0 = "Birney and Goldman"
exp_dna_out = "TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGC" \
          "TGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGTACGTACGTACGTACGTACGACTAT"
exp_split_dna = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGA","CATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTAC"]

exp_115 = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTC",  "ATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAG"]
exp_117 = ["ATAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTCG",  "TATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAGC"]


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
assert s5 == exp_117
print(s5)

def add_parity_trit(dna, ID):
    for i in range(len(dna)):
        dna[i] += b3_to_dna(generate_parity_trit(dna[i], i, ID), prev_char=dna[i][-1])
        assert len(dna[i]) == 115
    return dna

def polish_ends(dna):
    for i in range(len(dna)):
        dna[i] = add_final_trits(dna[i])
        assert len(dna[i]) == 117
    return dna

def encode(in_str, ID):
    ord_str = [ord(l) for l in in_str]
    b3_str = "".join(a2b3(ord_str))
    dna_str = build_initial_dna_string(b3_str)
    split_str = split_and_index(dna_str)
    out_str = polish_ends(add_parity_trit(split_str, ID))
    return out_str

print("xxxxxxxxxxxxxxxxx")
a = encode(s0, "12")
print(a)
assert a == exp_117