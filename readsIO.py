import gzip
from utils import rev_comp2, dna_map, dna_to_b3
from decode import b3_to_int, rev_comp_odd
from collections import Counter
# combine read pairs
# filter out pairs of incorrect length
# filter out read pairs with incorrect combos of start end char
# attempt to merge reads
file1 = "davos.join.fa.gz"
file2 = "davos_small.fa"
read_names = list()
read_seqs = list()

def load_reads(filename):
    with gzip.open(filename, 'rt') as f1:
        for line in f1.readlines():
            line = line.strip("\n")
            if not line.startswith(">"):
                read_seqs.append(line)
        return read_seqs
#reads = load_reads(file1)

def load_reads2(filename):
    with open(filename, 'rt') as f1:
        for line in f1.readlines(50000):
            read_seqs.append(line.strip("\n"))
        return read_seqs
reads = load_reads2(file2)
total_reads = len(reads)
print(total_reads)
filtered_reads = [read for read in reads if len(read) == 117]
filtered_reads_len = len(filtered_reads)
print(filtered_reads_len)
print(filtered_reads_len/total_reads)
filtered_reads2 = list()
bad = False
for read in filtered_reads:
    for i in range(len(read)-1):
        if read[i] == read[i+1]:
            bad = True
            break
    if not bad:
        filtered_reads2.append(read)
    bad = False
print("????")

print(len(filtered_reads2)/total_reads)
valid_starts = {'A', 'T'}
valid_ends = {'G', 'C'}
correct_ori_reads = list()
for read in filtered_reads2:
    if read[0] in valid_starts and read[-1] in valid_ends:
        correct_ori_reads.append(read[1:-1])
    else:
        read = rev_comp2(read)
        if read[0] in valid_starts and read[-1] in valid_ends:
            correct_ori_reads.append(read[1:-1])

def load_indexed_seqs(verified_dna):
    id_2 = [(dna_to_b3(seq[-15:], prev_char=seq[-16]), seq[:-15]) for seq in verified_dna]
    return dict(id_2)



correct_ori_reads_len = len(correct_ori_reads)
print(correct_ori_reads_len)
print(correct_ori_reads_len/total_reads)
indexed_b3 = load_indexed_seqs(correct_ori_reads)


def generate_parity_trit(i3, ID):
    i3 = "0" * (12 - len(i3)) + i3
    temp = ID + i3
    even_trits = [int(temp[x]) for x in range(len(temp)) if x % 2 == 0]
    parity_trit = sum(even_trits) % 3
    index = ID + i3 + str(parity_trit)
    return index

parity_counter = 0

def confirm_parity(index):

    file_id = index[:2]
    chunk_id = index[2:-1]
    if generate_parity_trit(chunk_id, file_id) == index:
        #logging.debug("Parity Trit OK!")
        out = (file_id, chunk_id)
    else:
        print("Parity error")
        #raise Exception("Parity Trit Error")
        out = ("9999", "9999")
    return out
def parse_index(indexed_b3):
    p = [(confirm_parity(x), indexed_b3[x]) for x in indexed_b3]
    print(dict(p)[("9999", "9999")])
    ### this step is killing duplicate chunk IDs!!!!1
    return dict(p)

def assign_to_files(parsed_index):
    files = dict()
    for key in parsed_index:
        file_id, chunk_id = key
        #print(file_id, chunk_id)
        if file_id not in files:
            files[file_id] = {chunk_id: list()}
            if file_id == "9999":
                print("if branch 1")
        else:
            if chunk_id not in files[file_id]:
                files[file_id][chunk_id] = list()
                if file_id == "9999":
                    print("else branch 2")
        files[file_id][chunk_id].append(parsed_index[key])
    # TODO: rewrite as generator
    print(files["9999"])
    return files


parsed_indexes = parse_index(indexed_b3)
file_id = Counter()
chunk_counter = Counter()
# for f_id, c_id in parsed_indexes:
#     file_id[f_id] += 1
#     if f_id == "00":
#         chunk_counter[c_id] += 1
print(file_id)
print(chunk_counter)
print(len(chunk_counter))
files = assign_to_files(parsed_indexes)



# file_00 = files["00"]
# chunk_counter = Counter()
# yu = [len(file_00[x]) for x in file_00]
# print("max=")
# print(max(yu))
# print(file_00[0])
# for k in file_00:
#     chunk_counter[k] += len(file_00[k])
# print(chunk_counter)
# print(min(chunk_counter.keys()), max(chunk_counter.keys()))
# missing = list()
# max_piece = 0
# for i in range(max(chunk_counter.keys())):
#     if i in chunk_counter:
#         max_piece = i
#     else:
#         break
# print(max_piece)
#
# print(max(chunk_counter.values()))
# print(len(file_00))
#
