import re
import string

from default_params import *
from pp.demux import main as demux_main
from utils import *


def preprocess(indir):
    r1_files = sorted(get_files(indir, r1_format))
    r2_files = sorted(get_files(indir, r2_format))
    ids = [re.findall(idpattern, os.path.basename(_))[0].strip(string.punctuation) for _ in r1_files]

    if demux_on:
        print("Demux is on, first start the demultiplexing process.")
        demux_dict['seqfile1'] = r1_files
        demux_dict['seqfile2'] = r2_files
        if not_overwrite_demux and os.path.isfile(demux_stats):
            print('demuxed files exist, pass it')
        else:
            r1_files, r2_files, ids, stats = demux_main(**demux_dict)
            stats_df = pd.DataFrame.from_dict(stats, orient='index')
            stats_df.loc[:, 'remaining reads/%'] = stats_df.loc[:, "remaining reads"] / stats_df.iloc[:, 1:4].sum(1)
            stats_df.to_csv(demux_stats, index=True)
        print("demultiplexing has been complete, indir will change to", demux_dir_samples)

    if not os.path.isfile(opath):
        write_manifest(opath=opath,
                       ids=ids,
                       r1_files=r1_files,
                       r2_files=r2_files, )
    print("manifest has been output to", opath)
    # 准备序列的输入
    raw_seq = import_data_with_manifest(opath)
    # 将序列信息导入qiime的环境,可另存为qza

    raw_seq_eval_vis = seq_eval(raw_seq,
                                n=n)
    raw_seq_eval_vis.save(os.path.join(odir,
                                       raw_seq_vis_path))

    join_params.update(qc_joined_params)
    join_params['raw_seq'] = raw_seq
    join_params['n'] = n
    joined_seq, joined_seq_eval_vis, \
    joined_qc_seq, joined_qc_eval_vis, joined_qc_stats = join_seqs(raw_seq,
                                                                   **join_params
                                                                   )
    os.makedirs(os.path.join(odir, 'preprocess'), exist_ok=True)
    raw_seq.save(os.path.join(odir, 'preprocess',
                              raw_seq_path))
    joined_qc_seq.save(os.path.join(odir, 'preprocess',
                                    joined_seq_path))
    joined_qc_eval_vis.save(os.path.join(odir, 'preprocess',
                                         joined_qc_seq_vis_path))
    joined_seq.save(os.path.join(odir, 'preprocess',
                                 joined_qc_seq_vis_path))
    joined_seq_eval_vis.save(os.path.join(odir, 'preprocess',
                                          joined_seq_vis_path))

    joined_qc_stats.view(pd.DataFrame).to_csv(os.path.join(odir, 'preprocess',
                                                           joined_qc_stats_tab_path), index=True)

    return raw_seq, joined_qc_seq