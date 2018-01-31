from utils import *
from functools import reduce
import logging

logging.basicConfig(level=logging.INFO)


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


def b3_to_int(b3):
    out = 0
    for i in range(len(b3)):
        rev = b3[::-1]
        out += (3 ** i) * int(rev[i])
    return out


def confirm_parity(index):
    file_id = index[:2]
    chunk_id = index[2:-1]
    if generate_parity_trit(chunk_id, file_id) == index:
        logging.debug("Parity Trit OK!")
    else:
        raise Exception("Parity Trit Error")
    out = (file_id, chunk_id)
    return out


def generate_parity_trit(i3, ID):
    i3 = "0" * (12 - len(i3)) + i3
    temp = ID + i3
    even_trits = [int(temp[x]) for x in range(len(temp)) if x % 2 == 0]
    parity_trit = sum(even_trits) % 3
    index = ID + i3 + str(parity_trit)
    return index


def merge_overlapping(dna1, dna2):
    if dna1[-75:] == dna2[:75]:
        if len(dna1) % 25000 == 0:
            logging.info("Merged {0} sequences".format(((len(dna1) - 100) / 25) + 1))
        return dna1 + dna2[75:]
    else:
        raise Exception("Could not merge dna strings \n dna1: {0} \n dna2: {1} ".format(dna1[-75:], dna2[:75]))


def strip_len_info(b3):
    str_len = b3_to_int(b3[-20:])
    b3 = b3[:-20]
    zeros_to_remove = len(b3) - str_len
    assert b3[-zeros_to_remove:] == "0" * zeros_to_remove
    b3 = b3[:-zeros_to_remove]
    return b3


def load_indexed_seqs(verified_dna):
    id_2 = [(dna_to_b3(seq[-15:], prev_char=seq[-16]), seq[:-15]) for seq in verified_dna]
    return dict(id_2)


def parse_index(indexed_b3):
    p = [(confirm_parity(x), indexed_b3[x]) for x in indexed_b3]
    return dict(p)

def rev_comp_odd(f1):
    for key in f1:
        if key % 2 != 0:
            f1[key] = rev_comp2(f1[key])
    return f1

def assign_to_files(parsed_index):
    files = dict()
    for key in parsed_index:
        file_id, chunk_id = key
        if file_id not in files:
            files[file_id] = {b3_to_int(chunk_id): parsed_index[key]}
        else:
            files[file_id][b3_to_int(chunk_id)] = parsed_index[key]
    # TODO: rewrite as generator
    return files


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


def merge_newest(f1):
    return reduce(merge_overlapping, map(lambda x: reduce(merge_overlapping, x), split_up(f1)))


def decode(in_filename, out_filename):
    with open(in_filename, 'r') as f:
        metamorph = [line.strip("\n") for line in f.readlines()]
    verified_dna = check_len_and_orientation(metamorph)
    indexed_b3 = load_indexed_seqs(verified_dna)
    parsed_index = parse_index(indexed_b3)
    f1 = assign_to_files(parsed_index)["12"]
    f1 = rev_comp_odd(f1)


    logging.info("Found {0} sequences, encoding approximately {1} characters".format(len(f1), len(f1) * 18))

    merged = merge_newest(f1)
    logging.info("Merge complete")
    logging.info("Begin dna to base3 conversion")
    b3_message = dna_to_b3(merged)

    ord_data = b3_to_ord(strip_len_info(b3_message))

    with open(out_filename, "wb") as f:
        f.write(ord_data)

#decode("out.jpg", "decoded.png")
