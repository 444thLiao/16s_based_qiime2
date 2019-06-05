import csv
import os
from glob import glob
from multiprocessing import Pool, cpu_count
from subprocess import check_call
import sys
import pandas as pd
try:
    from qiime2 import Artifact
    from qiime2.plugins.demux.visualizers import summarize as seq_summ_vis
    from qiime2.plugins.quality_filter.methods import q_score_joined
    from qiime2.plugins.vsearch.methods import join_pairs
except:
    pass
from skbio.io import read, write
from tqdm import tqdm


def data_parser(path, ft='csv', **kwargs):
    if type(path) != str and ft != 'metadatas':
        df = path.copy()
        if type(df) != pd.DataFrame:
            df = pd.DataFrame(df)
    else:
        if ft == 'csv':
            sniffer = csv.Sniffer()
            sniffer.preferred = [',', '\t', '|']
            dialect = sniffer.sniff(open(path, 'r').readline().strip('\n'))
            df = pd.read_csv(path, sep=dialect.delimiter, header=0, low_memory=False, **kwargs)
            # low_memory is important for sample name look like float and str. it may mistaken the sample name into some float.
        else:
            df = pd.read_excel(path, header=0, **kwargs)

    return df

def run_cmd(cmd, dry_run=False, log_file=None, **kwargs):
    outstream = None
    if type(log_file) == str:
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) != 0:
                outstream = open(log_file, 'a')
        if outstream is None:
            valid_path(log_file,check_ofile=True)
            outstream = open(log_file, 'w')
    elif log_file is None:
        outstream = sys.stdout
    else:
        outstream = log_file

    print(cmd, file=outstream)
    outstream.flush()
    if not dry_run:
        check_call(cmd,
                   shell=True,
                   executable="/usr/bin/zsh",
                   stdout=outstream,
                   stderr=outstream,
                   **kwargs)
        outstream.flush()
def valid_path(in_pth,
               check_size=False,
               check_dir=False,
               check_glob=False,
               check_odir=False,
               check_ofile=False):
    if type(in_pth) == str:
        in_pths = [in_pth]
    else:
        in_pths = in_pth[::]
    for in_pth in in_pths:
        in_pth = os.path.abspath(os.path.realpath(in_pth))
        if in_pth is None:
            continue
        if check_glob:
            query_list = glob(in_pth)
            if not query_list:
                raise Exception('Error because of input file pattern %s' % in_pth)
        if check_dir:
            if not os.path.isdir(in_pth):
                raise Exception("Error because %s doesn't exist" % in_pth)
        if check_size:
            if os.path.getsize(in_pth) <= 0:
                raise Exception("Error because %s does not contain content." % in_pth)
        if check_odir:
            if not os.path.isdir(in_pth):
                os.makedirs(in_pth, exist_ok=True)
        if check_ofile:
            odir_file = os.path.dirname(in_pth)
            if not os.path.isdir(odir_file):
                os.makedirs(odir_file, exist_ok=True)
    return True

def get_files(indir, p):
    f_pattern = p.strip('.')
    files = glob(os.path.join(indir, f_pattern), recursive=True)
    return list(sorted(files))


def write_manifest(opath, r1_files, r2_files, ids):
    valid_path(opath,check_ofile=True)

    template_text = "sample-id,absolute-filepath,direction\n"
    for r1, r2, sid in zip(r1_files, r2_files, ids):
        template_text += ','.join([sid, r1, 'forward']) + '\n'
        template_text += ','.join([sid, r2, 'reverse']) + '\n'
    with open(opath, 'w') as f1:
        f1.write(template_text)
    return opath


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
    try:
        raw_seq_eval_vis = seq_summ_vis(seq_data,
                                        n=n)
    except:
        import pdb;
        pdb.set_trace()
    # 可视化原始数据的质量评估,未joined,n为随机不放回采样的reads数量
    raw_seq_eval_vis = raw_seq_eval_vis[0]

    return raw_seq_eval_vis




def mv_seq(seq, opath, name_dict):
    seq = read(seq, format='fasta')
    with open(opath, 'w') as f1:
        for i in seq:
            pre_name = i.metadata['id']
            i.metadata['id'] = name_dict[pre_name]
            i.metadata['description'] = ''
            write(i, 'fasta', f1)


def save(obj, odir, name):
    obj.save(os.path.join(odir, name))


def parse_param(file, g):
    with open(file, 'r') as f1:
        exec(f1.read(), g)


def assign_work_pool(func, differ_args, num_thread):
    if num_thread == 0:
        # if number of thread equal to 0, then use all threads of computer.
        num_thread = cpu_count()
    with Pool(processes=num_thread) as pool:
        results = list(tqdm(pool.imap(func, differ_args), total=len(differ_args)))

    return results

# if __name__ == '__main__':
#     indir = '/home/liaoth/data2/16s/肾衰小鼠/raw_data'
#     opath = '/home/liaoth/data2/16s/qiime2_learn/shenshuai_manifest'
#     r1_format = '.*_1.fastq.gz'
#     idpattern = '(.*)_[12].fastq.gz'
#
#     write_manifest(indir,
#                    opath,
#                    r1_format,
#                    idpattern)
