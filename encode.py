import binascii
from random import randint
from utils import *

s0 = "Birney and Goldman"
exp_dna_out = "TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGC" \
          "TGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGTACGTACGTACGTACGTACGACTAT"
exp_split_dna = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGA","CATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTAC"]

exp_115 = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTC",  "ATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAG"]
exp_117 = ["ATAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTCG",  "TATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAGC"]




def build_initial_dna_string(b3_str):
    # calculate length for indexing
    str_to_encode = add_zeros(b3_str)

    return b3_to_dna(str_to_encode)


def split_and_index(dna):
    dna = chunk_str(dna)
    for i in range(len(dna)):
        if i % 2 == 1:
            dna[i] = rev_comp2(dna[i])
    return dna

def rand_base(bases):
    return bases[0]

def add_final_trits(dna):
    if dna.startswith("A"):
        dna = "T" + dna
    elif dna.startswith("T"):
        dna = "A" + dna
    else:
        dna = ["A", "T"][0] + dna #[randint(0, 10) % 2]
    if dna.endswith("C"):
        dna += "G"
    elif dna.endswith("G"):
        dna += "C"
    else:
        dna += ["C", "G"][0] #[randint(0, 10) % 2]
    #print(dna)
    return dna



def add_parity_trit(dna, ID):
    for i in range(len(dna)):
        dna[i] += b3_to_dna(generate_parity_trit(dna[i], i, ID), prev_char=dna[i][-1])
        assert len(dna[i]) == 115
    return dna

def polish_ends(dna):
    for i in range(len(dna)):
        assert len(dna[i]) == 115
        dna[i] = add_final_trits(dna[i])
        assert len(dna[i]) == 117, "Bad line: number {0} {1}".format(i+1, dna[i])
    return dna

def encode(in_str, ID):
    if not isinstance(in_str, bytes):
        ord_str = [ord(l) for l in in_str]
    else:
        ord_str = [int(x) for x in in_str]
    b3_str = "".join(a2b3(ord_str))
    dna_str = build_initial_dna_string(b3_str)
    split_str = split_and_index(dna_str)
    out_str = polish_ends(add_parity_trit(split_str, ID))
    return out_str

def create_dna_from_file(in_filename, out_filename):
    a = encode(s0, "12")
    print("xxxxxxxxxxxxxxxxx")
    assert a == exp_117
    with open(in_filename, 'rb') as f:
        metamorph = f.read()
        #m_ascii = metamorph.encode('ascii', 'ignore')
        #m_cleaned = m_ascii.decode('ascii')
    b = encode(metamorph, "12")
    with open(out_filename, 'w') as f2:
        f2.write("\n".join(b))

create_dna_from_file("background_img.png", "out.jpg")