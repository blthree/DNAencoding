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
            dna[i] = rev_comp(dna[i])
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

s0_ord = [ord(l) for l in s0]
s1 = "".join(a2b3(s0_ord))
s2 = build_initial_dna_string(s1)
assert s2 == exp_dna_out
#print(s2)
#print(exp_dna_out)
s5 = split_and_index(s2)
assert rev_comp("AATTGCA") == "TGCAATT"
ID = "12"
for i in range(len(s5)):
    s5[i] += b3_to_dna(generate_parity_trit(s5[i], i, ID), prev_char=s5[i][-1])
    assert len(s5[i]) == 115
#print(s5)
#print(exp_115)
assert s5 == exp_115
for i in range(len(s5)):
    s5[i] = add_final_trits(s5[i])
    assert len(s5[i]) == 117
# this is the final 117 bp dna string
# will be 125 bp with illumina adapters
assert s5 == exp_117
#print(s5)

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
    ord_str = [ord(l) for l in in_str]
    b3_str = "".join(a2b3(ord_str))
    dna_str = build_initial_dna_string(b3_str)
    split_str = split_and_index(dna_str)
    out_str = polish_ends(add_parity_trit(split_str, ID))
    return out_str


a = encode(s0, "12")
#print(a)
print("xxxxxxxxxxxxxxxxx")
assert a == exp_117
with open("pg5200.txt", 'r') as f:
    metamorph = "\n".join(f.readlines())
    m_ascii = metamorph.encode('ascii', 'ignore')
    m_cleaned = m_ascii.decode('ascii')

b = encode(m_cleaned, "12")
with open("out.txt", 'w') as f2:
    f2.write("\n".join(b))