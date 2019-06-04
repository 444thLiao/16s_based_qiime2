import os
from glob import glob
from os.path import dirname

from pp.demux import main as demux_main

basic_dir = dirname(__file__)
metadata = os.path.join(basic_dir, 'test', 'metadata.tab')
id_col = 'SampleID'
fb_col = 'Forward_Barcode'
rb_col = 'Reverse_Barcode'
fp_col = 'Forward_Primer'
rp_col = 'Reverse_Primer'
seqfile1, seqfile2 = sorted(glob(os.path.join(basic_dir,
                                              'test/seq_data2/test_seq*_1.fastq.gz')
                                 )), \
                     sorted(glob(os.path.join(basic_dir,
                                              'test/seq_data2/test_seq*_2.fastq.gz')))
# the order of glob output is random, be careful !!!!!!!!!!!!!!!!!!!!!!!!1
f1_files, f2_files, ids, stats = demux_main(metadata=metadata,
                                            id_col=id_col,
                                            rb_col=rb_col,
                                            fb_col=fb_col,
                                            rp_col=rp_col,
                                            fp_col=fp_col,
                                            seqfile1=seqfile1,
                                            seqfile2=seqfile2,
                                            output_dir_pre=os.path.join(basic_dir, 'test/seq2_demux/'),
                                            output_dir_samples=os.path.join(basic_dir, 'test/seq2_samples/'),
                                            attempt_read_orientation=True)
for f, s in sorted(stats.items(), key=lambda x: x[0]):
    print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.keys()))
    print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.values()))
