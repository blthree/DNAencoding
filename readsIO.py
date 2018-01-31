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
        for line in f1.readlines():
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
        #print("Parity error")
        #raise Exception("Parity Trit Error")
        out = ("9999", "9999")
    return out
def parse_index(indexed_b3):
    #p = [(confirm_parity(x), indexed_b3[x]) for x in indexed_b3]

    ### this step is killing duplicate chunk IDs!!!!1
    p = dict()
    for x in indexed_b3:
        file_id, chunk_id = confirm_parity(x)
        chunk_id = b3_to_int(chunk_id)
        seq = indexed_b3[x]
        if not file_id in p:
            p[file_id] = {chunk_id: list()}
        else:
            if chunk_id not in p[file_id]:
                p[file_id][chunk_id] = list()
        p[file_id][chunk_id].append(seq)

    return p


parsed_indexes = parse_index(indexed_b3)

file_id = Counter()

for f in parsed_indexes:
    for c in parsed_indexes[f]:
        file_id[f] += len(parsed_indexes[f][c])
print(file_id)
file_00 = parsed_indexes["00"]
chunk_counter = Counter()
for chunk in file_00:
    if len(file_00[chunk]) > 1:
        print(chunk, len(file_00[chunk]))
    chunk_counter[chunk] += len(file_00[chunk])
print(chunk_counter)

print(max(file_00.keys()))
for i in range(len(file_00)):
    if not i in file_00:
        print("Missing chunk # {0}".format(i))
