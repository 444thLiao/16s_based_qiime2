"""
For calculating the length of overlapping region.
not api,just demo...but too slow,abort.
"""
from subprocess import check_output
from glob import glob
from tqdm import tqdm
import os
from Bio import SeqIO
import gzip
from difflib import SequenceMatcher

from collections import defaultdict

length_dis = defaultdict(list)

for jf in tqdm(glob("/home/liaoth/data2/16s/shanghai_152/pipelines/joined_seq/*.gz")):
    sample_name = os.path.basename(jf).split('.f')[0]
    f1 = "/home/liaoth/data2/16s/shanghai_152/raw/%s.R1.fq.gz" % sample_name
    f2 = "/home/liaoth/data2/16s/shanghai_152/raw/%s.R2.fq.gz" % sample_name
    total_num_reads = int(check_output("zgrep -c '^+$' %s" % os.path.abspath(jf), executable='/usr/bin/zsh',
                                       shell=True))

    fq1 = SeqIO.parse(gzip.open(f1, 'rt'), format='fastq')
    fq2 = SeqIO.parse(gzip.open(f2, 'rt'), format='fastq')
    for rj in tqdm(SeqIO.parse(gzip.open(jf, 'rt'), format='fastq'), total=total_num_reads):

        n1 = n2 = n3 = n4 = 0
        while 1:
            r1 = next(fq1)
            r2 = next(fq2)

            if r1.description == rj.description or r2.description == rj.description:
                break

        m1 = SequenceMatcher(None, str(r1.seq), str(rj.seq))
        m2 = SequenceMatcher(None, str(r2.seq), str(rj.reverse_complement().seq))
        a, b = (m1.get_matching_blocks(),
                m2.get_matching_blocks())
        # if (a[0].size + b[0].size) < len(rj):
        #     print(a,b,len(rj))
        #     # m1 = SequenceMatcher(None, str(r1.seq), str(rj.reverse_complement().seq))
        #     # m2 = SequenceMatcher(None, str(r2.seq), str(rj.seq))
        #     # a, b = (m1.get_matching_blocks(),
        #     #         m2.get_matching_blocks())
        #     import pdb;pdb.set_trace()

        result = abs(len(rj) - a[0].size - b[0].size)
        length_dis[sample_name].append(result)

if __name__ == '__main__':
    pass
