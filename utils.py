import os
import re
from glob import glob
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize as seq_summ_vis
from qiime2.plugins.quality_filter.methods import q_score_joined
from qiime2.plugins.vsearch.methods import join_pairs


def write_manifest(indir, opath, r1_format, idpattern):
    template_text = "sample-id,absolute-filepath,direction\n"
    f_pattern = '*.gz'

    for fp in glob(os.path.join(indir, f_pattern), recursive=True):
        textid = re.findall(idpattern, os.path.basename(fp))[0]
        path = fp
        direction = 'forward' if re.match(r1_format, fp) else 'reverse'
        template_text += ','.join([textid, path, direction]) + '\n'

    with open(opath, 'w') as f1:
        f1.write(template_text)


def import_data_with_manifest(manifest):
    input_type = 'SampleData[PairedEndSequencesWithQuality]'
    input_format = 'PairedEndFastqManifestPhred64'
    try:
        raw_data = Artifact.import_data(input_type,
                                        manifest,
                                        view_type=input_format)
    except ValueError as e:
        if 'Decoded Phred score is out of range' in e.args[0]:
            input_format = 'PairedEndFastqManifestPhred33'
            raw_data = Artifact.import_data(input_type,
                                            manifest,
                                            view_type=input_format)
        else:
            return
    return raw_data


def load_classifier(c_path):
    # c_path = '/home/liaoth/data2/gg-13-8-99-nb-classifier.qza'
    if not c_path.endswith('.qza'):
        classifier = Artifact.load(c_path)
    else:
        pass
        # todo: implement a classifier training module.
    return classifier


def seq_eval(seq_data,
             n=10000):
    raw_seq_eval_vis = seq_summ_vis(seq_data,
                                    n=n)
    # 可视化原始数据的质量评估,未joined,n为随机不放回采样的reads数量
    raw_seq_eval_vis = raw_seq_eval_vis[0]

    return raw_seq_eval_vis


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
              ):
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
    joined_qc_seq, joined_qc_stats = q_score_joined(demux=joined_seq,
                                                    min_quality=min_quality,
                                                    quality_window=quality_window,
                                                    min_length_fraction=min_length_fraction,
                                                    max_ambiguous=max_ambiguous,
                                                    )
    joined_qc_eval_vis = seq_summ_vis(joined_qc_seq,
                                      n=n
                                      )
    joined_qc_df = joined_qc_stats.view(pd.DataFrame)


    return (joined_seq,
            joined_seq_eval_vis,
            joined_qc_seq,
            joined_qc_eval_vis,
            joined_qc_stats
            )


if __name__ == '__main__':
    indir = '/home/liaoth/data2/16s/肾衰小鼠/raw_data'
    opath = '/home/liaoth/data2/16s/qiime2_learn/shenshuai_manifest'
    r1_format = '.*_1.fastq.gz'
    idpattern = '(.*)_[12].fastq.gz'

    write_manifest(indir,
                   opath,
                   r1_format,
                   idpattern)
