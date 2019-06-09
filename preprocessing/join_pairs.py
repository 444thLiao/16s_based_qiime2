#################################################################################
####  用vsearch进行join PE sequencing data
####  20190609
####  by tianhua liao
#################################################################################
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import pandas as pd
from pp import constant, run_cmd, fileparser, valid_path
import multiprocessing as mp
import click
from glob import glob
from tqdm import tqdm


def exec_fun(args):
    func, args = args
    return func(*args)


join_params = constant.join_params


def join_with_vsearch(fwd_fp,
                      rev_fp,
                      fastq_out,
                      log_file,
                      truncqual: int = join_params['truncqual'],
                      minlen: int = join_params['minlen'],
                      maxns: int = join_params['maxns'],
                      allowmergestagger: bool = join_params['allowmergestagger'],
                      minovlen: int = join_params['minovlen'],
                      maxdiffs: int = join_params['maxdiffs'],
                      minmergelen: int = join_params['minmergelen'],
                      maxmergelen: int = join_params['maxmergelen'],
                      maxee: float = join_params['maxee'],
                      qmin: int = join_params['qmin'],
                      qminout: int = join_params['qminout'],
                      qmax: int = join_params['qmax'],
                      qmaxout: int = join_params['qmaxout']):
    cmd = ' '.join([constant.vsearch_path if constant.vsearch_path else '/usr/bin/vsearch',
                    '--fastq_mergepairs', fwd_fp,
                    '--reverse', rev_fp,
                    '--fastqout', fastq_out,
                    '--fastq_ascii', '33',  # todo: use 33 as default
                    '--fastq_minlen', str(minlen),
                    '--fastq_minovlen', str(minovlen),
                    '--fastq_maxdiffs', str(maxdiffs),
                    '--fastq_qmin', str(qmin),
                    '--fastq_qminout', str(qminout),
                    '--fastq_qmax', str(qmax),
                    '--fastq_qmaxout', str(qmaxout),
                    ])
    if truncqual is not None:
        cmd += ' --fastq_truncqual ' + str(truncqual)
    if maxns is not None:
        cmd += ' --fastq_maxns ' + str(maxns)
    if minmergelen is not None:
        cmd += ' --fastq_minmergelen ' + str(minmergelen)
    if maxmergelen is not None:
        cmd += ' --fastq_maxmergelen ' + str(maxmergelen)
    if maxee is not None:
        cmd += ' --fastq_maxee ' + str(maxee)
    if allowmergestagger:
        cmd += ' --fastq_allowmergestagger'
    sample_id = os.path.basename(fastq_out).replace('.fastq', '')

    run_cmd(cmd, dry_run=False, log_file=log_file)
    ori_num = os.popen("grep -c '^+$' %s " % fwd_fp).read().strip('\n')
    after_num = os.popen("grep -c '^+$' %s " % fastq_out).read().strip('\n')
    run_cmd("gzip %s" % fastq_out, dry_run=False, log_file=log_file)
    return ori_num, after_num, sample_id


def get_paths(files_input):
    if type(files_input) == str:
        return list(sorted(glob(files_input)))
    else:
        return files_input


@click.command()
@click.option("-m", '--metadata', "metadata", help='path of metadata')
@click.option("-o", '--output-dir', "output_dir", help='The directory you want to output to')
@click.option("-p", '--num_thread', "n_thread", default=0, show_default=True,
              help='The number of thread you want to use. 0 mean all threads')
@click.option("-r1", '--r1-files', "r1_file", type=str,
              help='forward reads you want to demuplexed. If you pass wildcard, it should use quote to quote it.')
@click.option("-r2", '--r2-files', "r2_file", type=str,
              help='Reversed reads you want to demuplexed. If you pass wildcard, it should use quote to quote it.')
@click.option("-id", '--sampleid', "sample_id", type=str,
              help='Sample id')
def main(metadata=None,
         output_dir=None,
         n_thread=None,
         r1_file=None,
         r2_file=None,
         sample_id=None,
         ):
    odir = os.path.join(output_dir,
                        constant.joined_seq_path)
    if os.path.isdir(output_dir):
        os.system("rm -fr %s" % odir)
    valid_path(odir, check_odir=1)
    log_file = os.path.join(odir, "joined_process.log")
    stats_file = os.path.join(odir, "joined_stats.csv")
    if n_thread == 0 or n_thread == -1:
        n_thread = mp.cpu_count()

    if metadata is None:
        data = pd.DataFrame([[r1_file,
                              r2_file]],
                            index=[sample_id],
                            columns=["path_R1", "path_R2"])

    else:
        data = fileparser(metadata)
        data = data.df
    all_args = []

    for sampleid, (fwd_fp, rev_fp) in data.iterrows():
        fastq_out = os.path.join(output_dir,
                                 constant.joined_seq_path,
                                 str(sampleid) + '.fastq')
        all_args.append((join_with_vsearch,
                         (fwd_fp,
                          rev_fp,
                          fastq_out,
                          log_file)
                         ))
    summary_df = pd.DataFrame(index=data.index, columns=["raw", "joined", "joined_ratio/%"])
    with mp.Pool(processes=n_thread) as thread_pool:
        for ori_num, aft_num, sid in tqdm(thread_pool.imap(exec_fun, all_args),
                                          total=len(all_args)):
            ratio = float(aft_num) / float(ori_num) * 100
            summary_df.loc[sid, :] = ori_num, aft_num, ratio
    summary_df.to_csv(stats_file, index=1, sep=',')


if __name__ == '__main__':
    main()

    # python3 preprocessing/join_pairs.py -m 'test/seq2_demux_data.csv' -p 5 -o 'test/seq2_output/'
    # python3 preprocessing/join_pairs.py -r1 test/seq2_demux/G445302653_R1.fastq -r2 test/seq2_demux/G445302653_R2.fastq -id G445302653 -o test/seq2_output/
