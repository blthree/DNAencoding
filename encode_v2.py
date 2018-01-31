from encode import *
keystream = {0: "002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221",
            1: "202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022",
            2: "221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100",
            3: "100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020"}

test_str = "Birney and Goldman"
exp_dna0 = "CGTGTACGACTCGACAGAGATGCAGAGACTAGACACTCTCTATACGCATCTCAGACACTCGTAGTCTACTGACATCGCACGAGAGAGAGCAGCGCACTCA"
def add_zeros2(b3_str, mod_num=25):
    b3_str_len = len(b3_str)
    b3_len_enc = int_to_b3(b3_str_len)
    b3_len_enc = "0" * (25 - len(b3_len_enc)) + b3_len_enc
    print(b3_len_enc)
    zeros_to_add = ((((b3_str_len + len(b3_len_enc)) // mod_num) + 1) * mod_num - (b3_str_len + len(b3_len_enc)))

    str_to_encode = b3_len_enc + b3_str + "0" * zeros_to_add
    assert (len(str_to_encode) % mod_num) == 0
    return str_to_encode
def build_initial_dna_string2(b3_str):
    # calculate length for indexing
    str_to_encode = add_zeros2(b3_str)

    return b3_to_dna(str_to_encode)

def encode(in_str, ID):
    if not isinstance(in_str, bytes):
        ord_str = [ord(l) for l in in_str]
    else:
        ord_str = [int(x) for x in in_str]
    b3_str = "".join(a2b3(ord_str))
    print(b3_str)
    dna_str = build_initial_dna_string2(b3_str)
    #split_str = split_and_index(dna_str)
    #out_str = polish_ends(add_parity_trit(split_str, ID))
    return dna_str
dna_str = encode(test_str, "12")
dna0 = dna_str[:100]
dna2base4 = str.maketrans({'A':'0', 'C':'1', 'G':'2', 'T':'3'})
base4_0 = dna0.translate(dna2base4)
key0 = [(int(base4_0[i+1]) - int(base4_0[i]))%3 for i in range(len(base4_0)-1)]
print(len(key0))
print(key0)
print(dna_str)
for i in range(len(key0)):
    key0[i] += int(keystream[0])
    if key0[i] > 3:
        key0[i] = key0[i] % 3
key1 = [(key0[i] + 1) for i in range(len(key0)-1)]
encrypted = list()

for i in range(len(base4_0)-1):
    if i == 0:
        encrypted.append(base4_0[i])
    else:
        new = (int(base4_0[i]) + key1[i-1]) % 4
        encrypted.append(new)
base4_1 = "".join([str(x) for x in encrypted])
base4_2dna = str.maketrans({'0':'A', '1':'C', '2':'G', '3':'T'})
dna2 = base4_1.translate(base4_2dna)
print(dna2)


