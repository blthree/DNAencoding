
dna_in = \
    ['ATAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGATACGTACGTACGTCG',
    'TATAGTCGTACGTACGTACGTACGTACGTACGTACTGTACAGAGTCACTCGTCATCGATACTCACAGCATGCTGCGTAGCAGCGTATCTCGCTGCGAGATGATACGTACGTACGAGC']
exp_dna = "TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTACGTACGTACGTACGTACGTACGACTAT"

exp_split_dna = ["TAGTATATCGACTAGTACAGCGTAGCATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGA", "CATCTCGCAGCGAGATACGCTGCTACGCAGCATGCTGTGAGTATCGATGACGAGTGACTCTGTACAGTACGTAC"]

exp_b3 = "20100202101010100021200012221110221201112000212210002212222212021100210122100110210111200021000000000000000000000000000010102"
huffman_decode = {'22201': 0, '22200': 85, '22122': 170, '22121': 127, '22120': 253, '22112': 52, '22111': 138, '22110': 41, '22102': 86, '22101': 42, '22100': 100, '22022': 44, '22020': 250, '22021': 132, '22012': 161, '22010': 98, '22002': 8, '22011': 34, '22001': 10, '22000': 149, '21222': 87, '21221': 21, '21220': 74, '21212': 36, '21210': 69, '21202': 177, '21211': 20, '21200': 213, '21201': 163, '21121': 229, '21122': 255, '21120': 197, '21112': 133, '21110': 252, '21111': 26, '21101': 173, '21102': 151, '21100': 82, '21022': 75, '21021': 37, '21011': 166, '21020': 191, '21012': 88, '21010': 63, '21001': 68, '21002': 150, '21000': 76, '20222': 4, '20221': 154, '20212': 234, '20220': 22, '20211': 162, '20210': 105, '20202': 102, '20201': 171, '20200': 104, '20122': 169, '20121': 196, '20120': 208, '20112': 84, '20111': 130, '20102': 146, '20110': 72, '20101': 16, '20100': 66, '20022': 24, '20012': 106, '20020': 223, '20021': 58, '20011': 137, '20010': 73, '20001': 101, '20002': 168, '12221': 181, '12222': 175, '20000': 251, '12220': 40, '12212': 140, '12211': 17, '12210': 83, '12202': 254, '12201': 240, '12200': 214, '12122': 53, '12112': 202, '12121': 25, '12120': 18, '12111': 247, '12110': 174, '12102': 112, '12101': 89, '12100': 210, '12012': 217, '12020': 248, '12021': 194, '12022': 182, '12011': 80, '12002': 79, '12010': 195, '12001': 12, '12000': 209, '11222': 165, '11221': 245, '11220': 2, '11212': 81, '11211': 38, '11202': 141, '11210': 211, '11200': 239, '11201': 95, '11122': 43, '11121': 224, '11112': 203, '11120': 145, '11110': 147, '11111': 19, '11101': 50, '11102': 136, '11100': 107, '11022': 134, '11021': 109, '11020': 153, '11002': 148, '11010': 205, '11011': 212, '11012': 54, '11000': 241, '11001': 156, '10222': 115, '10221': 116, '10220': 78, '10211': 67, '10212': 70, '10210': 178, '10202': 159, '10201': 142, '10200': 92, '10122': 48, '10120': 90, '10121': 218, '10112': 126, '10111': 39, '10102': 219, '10110': 167, '10101': 114, '10022': 172, '10100': 14, '10020': 120, '10021': 139, '10012': 160, '10011': 33, '10010': 179, '10002': 117, '10001': 225, '10000': 129, '02222': 183, '02220': 230, '02221': 35, '02210': 93, '02211': 6, '02212': 32, '02201': 56, '02202': 158, '02121': 185, '02122': 47, '02200': 143, '02111': 123, '02120': 204, '02112': 242, '02110': 111, '02102': 103, '02101': 108, '02100': 9, '02022': 65, '02020': 249, '2136': 2021, '02012': 180, '02001': 226, '02002': 144, '02010': 15, '02011': 57, '02000': 128, '01220': 135, '01221': 243, '01222': 190, '01212': 207, '01211': 77, '01210': 45, '01202': 91, '01201': 192, '01122': 186, '01200': 216, '01112': 97, '01120': 118, '01121': 246, '01111': 215, '01102': 51, '01110': 206, '01100': 184, '01101': 227, '01022': 233, '01021': 237, '01020': 188, '01012': 113, '01011': 49, '01010': 201, '01002': 155, '01000': 222, '01001': 231, '00222': 5, '00221': 27, '00212': 131, '00220': 164, '00211': 3, '00210': 46, '00201': 119, '00202': 28, '00200': 176, '00122': 23, '00121': 64, '00120': 157, '00112': 187, '00110': 244, '00111': 238, '00102': 96, '00101': 235, '00022': 60, '00100': 1, '00021': 110, '00011': 200, '00020': 221, '00012': 99, '00010': 31, '00002': 198, '00001': 193, '00000': 125, '222222': 124, '222221': 152, '222220': 122, '222212': 71, '222211': 94, '222210': 220, '222202': 29, '222201': 199, '222200': 61, '222122': 11, '222121': 228, '222120': 62, '222112': 55, '222111': 121, '222110': 7, '222102': 30, '222101': 232, '222100': 189, '222021': 59, '222022': 236}



def rev_comp(dna):
    complements = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    dna_out = ""
    for l in dna[::-1]:
        dna_out += complements[l]
    return dna_out


def check_len_and_orientation(dna):
    valid_starts = ['A', 'T']
    valid_ends = ['G', 'C']
    for i in range(len(dna)):
        seq = dna[i]
        if len(seq) == 117:
            if not seq[0] in valid_starts and seq[-1] in valid_ends:
                seq_rc = rev_comp(seq)
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
        print("Parity Trit OK!")
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

verified_dna = check_len_and_orientation(dna_in)
print(verified_dna)
indexed_dna = dict()
indexed_b3 = dict()
for seq in verified_dna:
    indexed_dna[seq[-15:]] = seq[:-15]
    indexed_b3[dna_to_b3(seq[-15:], prev_char=seq[-16])] = seq[:-15]
print(indexed_dna)
print(indexed_b3)
dna_map = {'A': {'C': '0', 'G': '1', 'T': '2'},
           'C': {'G': '0', 'T': '1', 'A': '2'},
           'G': {'T': '0', 'A': '1', 'C': '2'},
           'T': {'A': '0', 'C': '1', 'G': '2'}}

def merge_overlapping(dna1, dna2):
    if dna1[-75:] == dna2[:75]:
        return dna1 + dna2[-25:]
    else:
        raise Exception("Could not merge dna strings")

parsed_index = {}
for x in indexed_b3:
    parsed_index[confirm_parity(x)] = indexed_b3[x]

#print(parsed_index)
files = dict()
for key in parsed_index:
    file_id, chunk_id = key
    if file_id not in files:
        files[file_id] = {int(chunk_id): parsed_index[key]}
    else:
        files[file_id][int(chunk_id)] = parsed_index[key]
#print(files)# good as of here
f1 = files["12"]
for key in f1:
    if key % 2 != 0:
        f1[key] = rev_comp(f1[key])
#print(f1)
print(f1[0][-75:])
print(f1[1][:75])
merged = merge_overlapping(f1[0], f1[1])
print(merged)
assert merged == exp_dna
b3_message = dna_to_b3(merged)
print(b3_message)
assert b3_message == exp_b3

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
ord_data = b3_to_ord(strip_len_info(b3_message))
print(ord_data)
#print(strip_len_info(b3_message))
decoded = [chr(o) for o in ord_data]
print("".join(decoded))


