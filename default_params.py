# 输入
import os

_p = os.path.abspath(__file__)
_p_dir = os.path.dirname(_p)
indir = os.path.join(_p_dir, 'test', 'seq_data')
r1_format = '.*_1.fastq.gz'
r2_format = '.*_2.fastq.gz'
idpattern = '(.*).[12].fastq.gz'

# 输出
odir = '/tmp/test_result'
opath = '/tmp/test_result/seq_manifest'
log_name = "cmd_log.txt"
############################################################
# dir
preprocess_path = "preprocess"
############################################################
# 输出文件名集合
raw_seq_path = 'raw_data'
joined_seq_path = 'joined_seq'
joined_qc_seq_path = 'joined_qc_seq'

summarized_stats_path = "summary_stats.csv"
raw_seq_vis_path = 'unjoined_seq_eval'
joined_seq_vis_path = 'joined_seq_eval'
joined_qc_seq_vis_path = 'joined_qc_seq_eval'

profiled_tab_path = 'profiling.qza'
representative_sequence_path = 'rep_fa.qza'
process_stats_path = 'profiling_stats.qza'

root_tree_path = 'rep_rooted_tree.tab'
tax_tab = 'rep_sintax.tab'

demux_dir_pre = os.path.join(odir, 'demux', 'pre')
demux_dir_samples = os.path.join(odir, 'demux', 'samples')
demux_stats = os.path.join(demux_dir_pre, '..', 'stats.csv')
############################################################
# demux部分
demux_on = False
not_overwrite_demux = True
demux_dict = dict(
    metadata=os.path.join(_p_dir, 'test', 'metadata.tab'),
    id_col='SampleID',
    fb_col='Forward_Barcode',
    rb_col='Reverse_Barcode',
    fp_col='Forward_Primer',
    rp_col='Reverse_Primer',
    attempt_read_orientation=True,
    output_dir_pre=demux_dir_pre,
    output_dir_samples=demux_dir_samples,
    num_thread=0  # all threads
)

# 序列评估 可视化部分参数
n = 10000
# join 部分参数
join_params = dict(
    truncqual=None,
    minlen=1,
    maxns=None,
    allowmergestagger=True,
    minovlen=10,
    maxdiffs=10,
    minmergelen=None,
    maxmergelen=None,
    maxee=None,
    qmin=0,
    qminout=0,
    qmax=41,
    qmaxout=41
)
# join 序列评估
qc_joined_params = dict(
    min_quality=4,  # default
    quality_window=3,  # default
    min_length_fraction=0.75,  # default
    max_ambiguous=0,  # default
)

deblur_args = dict(
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
)

# pipeliens args
dada2_args = dict(
    # dada2
    trunc_len_f=140,
    trunc_len_r=150,
    n_threads=0,  # all threads
    trunc_q=2,  # default
    n_reads_learn=1000000,  # default
    max_ee=2.0,  # default

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

############################################################
#software path
vsearch_path = '/home/liaoth/tools/vsearch-2.6.2/bin/vsearch'
