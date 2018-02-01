from utils import dna_to_b3, huffman_decode, b3_to_ord
from decode import b3_to_int
dna0 = 'CGTGTACGACTCGACAGTGACACTCGATCAGTACATCGTCTGTCTGCGCTGTCTGTGTCTGTGTCACTACTATGATATGTATCGCGATAGAGCGCGCTCG'
keystream = {0: "0002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221",
            1: "202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022",
            2: "221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100",
            3: "100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020"}


dna2base4 = str.maketrans({'A':'0', 'C':'1', 'G':'2', 'T':'3'})
base4_2_dna = str.maketrans({'0':'A', '1':'C', '2':'G', '3':'T'})
print(dna0)
base4_0 = [int(x) for x in dna0.translate(dna2base4)]
kstrits = [int(x) for x in keystream[0]]
L = len(kstrits)
# this code adapted from the author's perl script, it did not translate well, but it works
# need to change keystream based on chunk_id tho
for i in range(1, L):
    j = L-i
    base4_0[j] = (base4_0[j] - base4_0[j-1]) % 4
print("".join(str(x) for x in base4_0))
for i in range(1, L):
    base4_0[i] = (base4_0[i]-1) % 3
print("".join(str(x) for x in base4_0))

for i in range(1, L):
    # add to encrypt, subtract to decrypt
    base4_0[i] = (base4_0[i] - kstrits[i]) % 3
print("".join(str(x) for x in base4_0))
for i in range(1, L):
    base4_0[i] = (base4_0[i] +1) % 4
print("".join(str(x) for x in base4_0))
for i in range(1, L-1):
    base4_0[i] = (base4_0[i] + base4_0[i-1]) % 4
print("".join(str(x) for x in base4_0))
temp6 = "".join([str(x) for x in base4_0])
temp6 = temp6.translate(base4_2_dna)
print(temp6)

b3_out = dna_to_b3(temp6)
print(b3_out)
len_info = b3_to_int(b3_out[:25])
print(len_info)
start = b3_to_ord(b3_out[25:])
deco = [chr(x) for x in start]
print("".join(deco))
print(len_info/25)