from pipelines import selective_p, tax_assign_qiime2, g_tree
from utils import *


def preprocess():
    write_manifest(indir,
                   opath,
                   r1_format,
                   idpattern)
    # 准备序列的输入
    raw_seq = import_data_with_manifest(opath)
    # 将序列信息导入qiime的环境,可另存为qza

    raw_seq_eval_vis = seq_eval(raw_seq,
                                n=n)
    raw_seq_eval_vis.save(os.path.join(odir,
                                       raw_seq_vis_path))

    joined_seq, joined_seq_eval_vis, \
    joined_qc_seq, joined_qc_eval_vis, joined_qc_stats = join_seqs(raw_seq,
                                                                   minlen=minlen,
                                                                   allowmergestagger=allowmergestagger,
                                                                   minovlen=minovlen,  # default
                                                                   maxdiffs=maxdiffs,  # default
                                                                   n=n,
                                                                   min_quality=min_quality,
                                                                   quality_window=quality_window,
                                                                   min_length_fraction=min_length_fraction,
                                                                   max_ambiguous=max_ambiguous,
                                                                   )
    raw_seq.save(os.path.join(odir, 'raw_data'))
    joined_qc_seq.save(os.path.join(odir,
                                    joined_seq_vis_path))
    joined_qc_eval_vis.save(os.path.join(odir,
                                         joined_qc_seq_vis_path))
    joined_seq.save(os.path.join(odir,
                                 joined_qc_seq_vis_path))
    joined_seq_eval_vis.save(os.path.join(odir,
                                          joined_seq_vis_path))

    joined_qc_stats.view(pd.DataFrame).to_csv(os.path.join(odir,
                                                           joined_qc_stats_tab_path), index=True)

    return raw_seq, joined_qc_seq


def run_pipelines(p, pipelines_args):
    pre_ = 'sOTU'
    tab, rep, tab_vis, seq_vis, stats_df = selective_p(p, pipelines_args)

    p_tab = tab.view(pd.DataFrame)
    name_dict = dict(zip(p_tab.columns, [pre_ + str(_) for _ in range(p_tab.shape[1])]))
    p_tab.columns = [name_dict[_] for _ in p_tab.columns]

    seq_f = glob(str(rep._archiver.data_dir) + '/*')
    if len(seq_f) != 1:
        raise Exception
    seq_f = seq_f[0]

    # output part
    tab.save(os.path.join(odir,
                              profiled_tab_path.format(prefix=p)))
    p_tab.to_csv(os.path.join(odir,
                              profiled_tab_path.format(prefix=p)),
                 sep='\t' if profiled_tab_path.endswith('.tab') else ',',
                 index=True)
    mv_seq(seq_f,
           opath=os.path.join(odir,
                              representative_sequence_path.format(prefix=p)),
           name_dict=name_dict)
    rep.save(os.path.join(odir,
                              representative_sequence_path.format(prefix=p)))
    tab_vis.save(os.path.join(odir,
                              profiled_tab_vis_path.format(prefix=p)))
    seq_vis.save(os.path.join(odir,
                              representative_sequence_vis_path.format(prefix=p)))
    stats_df.to_csv(os.path.join(odir,
                                 process_stats_path.format(prefix=p)))

    return tab, p_tab, rep


def after_otu(args):
    # assigning taxonomy, perform alpha,beta diversity
    tax_tab = tax_assign_qiime2(**args)
    rooted_tree = g_tree(**args)


if __name__ == '__main__':
    # 输入
    indir = '/home/liaoth/data2/16s/qiime2_learn/test_data'
    r1_format = '.*_1.fastq.gz'
    idpattern = '(.*)_[12].fastq.gz'
    metadata = '/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/metadata.tsv'
    # metadata = '/home/liaoth/data2/16s/qiime2_learn/meta_11.16.tsv'

    # 输出
    odir = '/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/'
    opath = os.path.join(odir, 'seq_manifest')
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
    ## 跑命令
    raw_seq, joined_qc_seq = preprocess()
    print("预处理完成,完成原始序列评估 与 joined, 去污染,去chimera,fix orientation")
    pipelines_args['deblur_input'] = joined_qc_seq
    pipelines_args['dada2_input'] = raw_seq

    p = 'dada2'
    tab, p_tab, rep = run_pipelines(p, pipelines_args)
    # after_otu_args['rep'] = rep

    # after_otu(after_otu_args)
