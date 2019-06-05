from qiime2.plugins.demux.visualizers import summarize as seq_summ_vis
from qiime2.plugins.quality_filter.methods import q_score_joined
from qiime2.plugins.vsearch.methods import join_pairs

def join_seqs(raw_data,
              minlen=100,
              allowmergestagger=True,
              minovlen=10,
              maxdiffs=10,
              n=10000,
              min_quality=4,  # default
              quality_window=3,  # default
              min_length_fraction=0.75,  # default
              max_ambiguous=0,  # default
              **kwargs
              ):
    print("Join_pairs starting.......")
    joined_seq = join_pairs(demultiplexed_seqs=raw_data,
                            minlen=minlen,
                            allowmergestagger=allowmergestagger,
                            minovlen=minovlen,  # default
                            maxdiffs=maxdiffs,  # default
                            )
    joined_seq = joined_seq[0]
    joined_seq_eval_vis = seq_summ_vis(joined_seq,
                                       n=n
                                       )
    print("Quality control at joined seq starting.......")
    joined_qc_seq, joined_qc_stats = q_score_joined(demux=joined_seq,
                                                    min_quality=min_quality,
                                                    quality_window=quality_window,
                                                    min_length_fraction=min_length_fraction,
                                                    max_ambiguous=max_ambiguous,
                                                    )
    joined_qc_eval_vis = seq_summ_vis(joined_qc_seq,
                                      n=n
                                      )

    return (joined_seq,
            joined_seq_eval_vis[0],
            joined_qc_seq,
            joined_qc_eval_vis[0],
            joined_qc_stats
            )