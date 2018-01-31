import gzip
from utils import rev_comp2, dna_map, dna_to_b3
from collections import Counter
# combine read pairs
# filter out pairs of incorrect length
# filter out read pairs with incorrect combos of start end char
# attempt to merge reads
file1 = "davos.join.fa.gz"

read_names = list()
read_seqs = list()

def load_reads(filename):
    with gzip.open(filename, 'rt') as f1:
        for line in f1.readlines():
            line = line.strip("\n")
            if line.startswith(">"):
                read_names.append(line)
            else:
                read_seqs.append(line)
        return dict(zip(read_names, read_seqs))
reads = load_reads(file1)

total_reads = len(reads)
print(total_reads)
#read_lens = [len(x) for x in reads.values()]
#print(sum(read_lens)/ len(read_lens))
filtered_reads = [reads[key] for key in reads if len(reads[key]) == 117]
filtered_reads_len = len(filtered_reads)
print(filtered_reads_len)
print(filtered_reads_len/total_reads)
filtered_reads2 = list()
bad = False
for read in filtered_reads:
    for i in range(len(read)-1):
        if read[i] == read[i+1]:
            bad = True
    if bad:
        filtered_reads2.append(read)
    bad = False


print(len(filtered_reads2)/total_reads)
valid_starts = {'A', 'T'}
valid_ends = {'G', 'C'}
correct_ori_reads = list()
for read in filtered_reads2:
    if read[0] in valid_starts and read[-1] in valid_ends:
        correct_ori_reads.append(read[1:-2])
    else:
        read = rev_comp2(read)
        if read[0] in valid_starts and read[-1] in valid_ends:
            correct_ori_reads.append(read[1:-2])
correct_ori_reads_len = len(correct_ori_reads)
print(correct_ori_reads_len)
print(correct_ori_reads_len/total_reads)
b3_reads = [dna_to_b3(read) for read in correct_ori_reads]
whole_indexes = [read[-15:] for read in b3_reads]

def generate_parity_trit(i3, ID):
    i3 = "0" * (12 - len(i3)) + i3
    temp = ID + i3
    even_trits = [int(temp[x]) for x in range(len(temp)) if x % 2 == 0]
    parity_trit = sum(even_trits) % 3
    index = ID + i3 + str(parity_trit)
    return index
def confirm_parity(index):
    file_id = index[:2]
    chunk_id = index[2:-1]
    if generate_parity_trit(chunk_id, file_id) == index:
        #logging.debug("Parity Trit OK!")
        out = (file_id, chunk_id)
    else:
        #raise Exception("Parity Trit Error")
        out = ("9999", "9999")

    return out

parsed_indexes = [confirm_parity(x) for x in whole_indexes]
file_id = Counter()
for f_id, c_id in parsed_indexes:
    file_id[f_id] += 1
print(file_id)

