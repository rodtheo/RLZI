import pandas as pd
import sys
from Bio import SeqIO
import math

filter_threshold = 50

out_chr_info = [["Chr", "Start", "End", "fill", "species", "size", "color"]]
offsets = [0]
offset = 0
mapchr = {}
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    out_chr_info.append([rec.id, 1, len(rec.seq), "969696","1", "12", "252525"])
    offset += len(rec.seq) + 1
    offsets.append(offset)
    offset += len(rec.seq) + 1
    offsets.append(offset)
for rec in SeqIO.parse(sys.argv[2], 'fasta'):
    out_chr_info.append([rec.id, 1, len(rec.seq), "969696", "2", "12", "252525"])

with open(sys.argv[4], "w") as outfchr:
    for r in out_chr_info:
        outfchr.write('{}\n'.format('\t'.join([str(x) for x in r])))


df = pd.read_csv(sys.argv[3], sep="\t", header=None)

# offset = offset + 1

header = "Species_1\tStart_1\tEnd_1\tSpecies_2\tStart_2\tEnd_2\tfill"
out_all = [header]

i = 1
for index, row in df.iterrows():
    _, pos_ref, len_phrase, _, chr = row
    if int(len_phrase) < -1:
        i += len_phrase
        continue
    if ((chr % 2) == 0):
        pos_ref_s = pos_ref - int(offsets[int(chr)-1])
        pos_ref_e = pos_ref + len_phrase - 1 - int(offsets[int(chr)-1])

        pos_ref_s_ok = int(offsets[int(chr)]) - (pos_ref + len_phrase - 1)
        pos_ref_e_ok = int(offsets[int(chr)]) - pos_ref
        chr = int(math.ceil(chr/2))
        out = [chr, pos_ref_s_ok + 1, pos_ref_e_ok , '1', i, i+len_phrase, "FF5733"]
    else:
        pos_ref_s = pos_ref - int(offsets[int(chr)-1])
        chr = int(math.ceil(chr/2)+1)
        out = [chr, pos_ref_s, pos_ref_s + len_phrase , '1', i, i+len_phrase, "cccccc"]
    i += len_phrase - 1
    # print("\t".join([str(x) for x in out]))
    out_all.append("\t".join([str(x) for x in out]))
# print(i)

fout = open(sys.argv[6],'w')
for nl, l in enumerate(out_all):
    if nl == len(out_all):
        fout.write(l)
    else:
        fout.write(l+'\n')

header = "Species_1\tStart_1\tEnd_1\tSpecies_2\tStart_2\tEnd_2\tfill"
out_filtered = [header]
i = 1
for index, row in df.iterrows():
    _, pos_ref, len_phrase, _, chr = row
    if int(len_phrase) <= filter_threshold:
        i += len_phrase
        continue
    if ((chr % 2) == 0):
        pos_ref_s = pos_ref - int(offsets[int(chr)-1])
        pos_ref_e = pos_ref + len_phrase - 1 - int(offsets[int(chr)-1])

        pos_ref_s_ok = int(offsets[int(chr)]) - (pos_ref + len_phrase - 1)
        pos_ref_e_ok = int(offsets[int(chr)]) - pos_ref
        chr = int(math.ceil(chr/2))
        out = [chr, pos_ref_s_ok + 1, pos_ref_e_ok , '1', i, i+len_phrase, "FF5733"]
    else:
        pos_ref_s = pos_ref - int(offsets[int(chr)-1])
        chr = int(math.ceil(chr/2)+1)
        out = [chr, pos_ref_s, pos_ref_s + len_phrase , '1', i, i+len_phrase, "cccccc"]
    i += len_phrase - 1
    # print("\t".join([str(x) for x in out]))
    out_filtered.append("\t".join([str(x) for x in out]))

fout = open(sys.argv[5],'w')
for nl, l in enumerate(out_filtered):
    if nl == len(out_filtered):
        fout.write(l)
    else:
        fout.write(l+'\n')
