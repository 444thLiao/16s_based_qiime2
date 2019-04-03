from qiime2.plugins.dada2.methods import denoise_paired
from qiime2.plugins.demux.visualizers import summarize as seq_summ_vis
from qiime2.plugins.feature_table.visualizers import summarize as ft_vis
from qiime2.plugins.feature_table.visualizers import tabulate_seqs as ft_tabseq
from qiime2.plugins.metadata.visualizers import tabulate

from utils import *


def dada2_pipelines(raw_seq_data,
                    trunc_len_f=140,
                    trunc_len_r=150,
                    n_threads=0,  # all threads
                    trunc_q=2,  # default
                    n_reads_learn=1000000,  # default
                    max_ee=2.0  # default
                    ):

    tab, rep, stats = denoise_paired(raw_seq_data,
                                     trunc_len_f=trunc_len_f,
                                     trunc_len_r=trunc_len_r,
                                     n_threads=n_threads,
                                     trunc_q=trunc_q,
                                     n_reads_learn=n_reads_learn,
                                     max_ee=max_ee
                                     )
    # 进行dada2
    dada2_tab_vis = ft_vis(tab)
    dada2_seq_vis = ft_tabseq(rep)
    dada2_stats_vis = tabulate(stats.view(qiime2.Metadata))

    return (raw_seq_eval_vis,
            tab,
            rep,
            stats,
            dada2_tab_vis,
            dada2_seq_vis,
            dada2_stats_vis)
