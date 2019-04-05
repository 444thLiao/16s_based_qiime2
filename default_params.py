# 输入
import os

_p = os.path.abspath(__file__)
_p_dir = os.path.dirname(_p)
indir = os.path.join(_p_dir, 'test', 'seq_data')
r1_format = '.*.1.fastq.gz'
r2_format = '.*.2.fastq.gz'
idpattern = '(.*).[12].fastq.gz'
metadata = os.path.join(_p_dir, 'metadata.tsv')

# 输出
odir = '/tmp/test_result'
opath = '/tmp/test_result/seq_manifest'

############################################################
# 文件名集合
raw_seq_vis_path = 'unjoined_seq_eval'
joined_seq_vis_path = 'joined_seq_eval'
joined_qc_seq_vis_path = 'joined_qc_seq_eval'
joined_qc_stats_tab_path = 'joined_qc_stats.csv'

profiled_tab_path = '{prefix}_profiling.tab'
profiled_tab_vis_path = '{prefix}_profiling'
representative_sequence_path = '{prefix}_rep.fa'
representative_sequence_vis_path = '{prefix}_rep'
process_stats_path = '{prefix}_profiling_stats.csv'

root_tree_path = '{prefix}_rep_rooted_tree.tab'
tax_tab = '{prefix}_rep_sintax.tab'

############################################################
# demux部分
demux_on = False
demux_dict = dict(
    metadata_file=os.path.join(_p_dir, 'test', 'metadata.tab'),
    id_col='SampleID',
    fb_col='Forward_Barcode',
    rb_col='Reverse_Barcode',
    fp_col='Forward_Primer',
    rp_col='Reverse_Primer',

)
# 序列评估 可视化部分参数
n = 10000
# join 部分参数
minlen = 100
allowmergestagger = True
minovlen = 10
maxdiffs = 10
# join 序列评估
min_quality = 4  # default
quality_window = 3  # default
min_length_fraction = 0.75  # default
max_ambiguous = 0  # default

pipelines_args = dict(
    # dada2
    trunc_len_f=140,
    trunc_len_r=150,
    n_threads=0,  # all threads
    trunc_q=2,  # default
    n_reads_learn=1000000,  # default
    max_ee=2.0,  # default

    # deblur
    trim_length=250,
    sample_stats=True,
    mean_error=0.005,  # default
    indel_prob=0.01,  # default
    indel_max=3,  # default
    min_reads=10,  # default
    min_size=2,  # default
    jobs_to_start=7,
    hashed_feature_ids=True,  # default

    # vsearch-otu

)

# Not implement yet
## after profiling
# after_otu_args = dict(
#     # classify
#     classifier_pth='/home/liaoth/data2/gg-13-8-99-nb-classifier.qza',
#     n_jobs=-2,
#     confidence=-1,
#     read_orientation=None,
#     # mafft tree
#     n_threads=0,
#     mask_max_gap_frequency=1.0,
#     mask_min_conservation=0.4,
# )
