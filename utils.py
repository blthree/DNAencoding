def load_huffman_code():
    with open("View_huff3.cd.new.correct.txt", "rb") as f:
        byte_base3 = dict()
        for line in f.readlines():
            line = line.split(b'\t')
            letter = line[1]
            byte = int(line[2])
            base3 = bytes.decode(line[3])
            byte_base3[byte] = base3
        return byte_base3
print(load_huffman_code())