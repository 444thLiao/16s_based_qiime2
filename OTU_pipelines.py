from qiime2.plugins.deblur.methods import denoise_16S

def deblur_p(joined_qc_seq):
    deblur_tab, deblur_rep, deblur_stats = denoise_16S(demultiplexed_seqs=joined_qc_seq,
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