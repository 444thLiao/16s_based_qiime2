from utils import *
import pandas as pd

def preprocess():
    write_manifest(indir,
                   opath,
                   r1_format,
                   idpattern)
    # 准备序列的输入
    seq_raw_data = import_data_with_manifest(opath)
    # 将序列信息导入qiime的环境,可另存为qza

    raw_seq_eval_vis = seq_eval(seq_raw_data,
                                n=n)
    raw_seq_eval_vis.save(os.path.join(odir,
                                       raw_seq_eval_vis))

    joined_seq, joined_seq_eval_vis, \
    joined_qc_seq, joined_qc_eval_vis, joined_qc_stats = join_seqs(seq_raw_data,
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
    joined_seq_eval_vis.save(os.path.join(odir,
                                      joined_seq_vis_path))
    joined_qc_eval_vis.save(os.path.join(odir,
                                      joined_qc_seq_vis_path))
    joined_qc_stats.to_csv(os.path.join(odir,
                                      joined_qc_seq_vis_path),index=True)
def run_pipelines(p):
    if p not in ['dada2',
                 'deblur',
                 'otu']:
        return
    else:
        pass

def after_otu():
    pass


if __name__ == '__main__':
    # 输入
    indir = '/home/liaoth/data2/16s/肾衰小鼠/raw_data'
    r1_format = '.*_1.fastq.gz'
    idpattern = '(.*)_[12].fastq.gz'
    metadata = '/home/liaoth/data2/16s/qiime2_learn/meta_11.16.tsv'

    # 输出
    odir = '/home/liaoth/data2/16s/qiime2_learn/shenshuai/'
    opath = os.path.join(odir, 'seq_manifest')
    ############################################################
    # 文件名集合
    raw_seq_vis_path = 'unjoined_seq_eval_vis'
    joined_seq_vis_path = 'joined_seq_eval_vis'
    joined_qc_seq_vis_path = 'joined_qc_seq_eval_vis'
    joined_qc_stats_tab_path = 'joined_qc_stats.csv'


    # 序列评估 可视化部分参数
    n = 10000
    # join 部分参数
    minlen = 100
    allowmergestagger = True
    minovlen = 10
    maxdiffs = 10
    # join 序列评估


    min_quality = 4,  # default
    quality_window = 3,  # default
    min_length_fraction = 0.75,  # default
    max_ambiguous = 0,  # default
    ## 跑命令
    preprocess()
    print("预处理完成,完成原始序列评估 与 joined, 去污染,去chimera,fix orientation")
