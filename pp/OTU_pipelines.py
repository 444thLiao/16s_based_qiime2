import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from glob import glob
from tqdm import tqdm
import gzip
from Bio import SeqIO
from subprocess import check_call
from default_file_structures import *


def task_screen(indir, joined_seq_path):
    outdir = os.path.join(indir, screen_dir)
    os.makedirs(outdir, exist_ok=True)
    fqs = os.path.join(indir, joined_seq_path, '*.fastq.gz')
    cmd = "fastq_screen {input} --outdir {outdir} --nohits --aligner bowtie2".format(input=fqs,
                                                                                     outdir=outdir
                                                                                     )
    check_call(cmd, shell=True)


def task_merge(indir, dirpath):
    # dirpath could be 1. screen_dir or 2. joined_seq_path
    indir = os.path.join(indir, dirpath, '*.fastq.gz')
    odir = os.path.join(indir, merged_dir)
    os.makedirs(odir, exist_ok=True)
    opath = os.path.join(odir, merged_file)
    # splitted by sample id; merge sequencing file joined after join_pairs.py
    all_fqs = sorted(glob(os.path.join(indir, '*.gz')))
    with open(opath, 'w') as f1:
        count = 0
        for fq in tqdm(all_fqs):
            stream = SeqIO.parse(gzip.open(fq, 'rt'), format='fastq')
            sid = os.path.basename(fq).replace('fastq.gz', '')
            for read in stream:
                read.id = read.description = read.name = ''
                read.id = '{sid}_{num};barcodelabel={sid}'.format(sid=sid,
                                                                  num=str(count))
                count += 1
                SeqIO.write(read, f1, format='fastq')


def task_filter(indir):
    input = os.path.join(indir, merged_dir,merged_file)
    output = os.path.join(indir, merged_dir,filtered_file)
    cmd = "vsearch --fastx_filter {input} --fastq_maxee 1.0 --fastq_trunclen 240 --fastaout {output}".format(input=input,
                                                                                                             output=output)
    check_call(cmd, shell=True)


def task_derep(indir):
    input = os.path.join(indir, merged_dir, merged_file)
    output = os.path.join(indir, merged_dir, filtered_file)
    cmd1 = "vsearch --derep_fulllength {input} --output {output} -sizeout".format(input=input,output=output)
    cmd2 = "vsearch --sortbysize  {input} --output {output} --minsize 1".format(input=input,output=output)
    check_call(cmd1)
    check_call(cmd2)

def task_remove_chimer():
    cmd = """vsearch --cluster_size splited/derep.fa \
	--id 0.97 --sizeout --fasta_width 0 \
	--uc v_analysis/all.preclustered.uc \
	--centroids v_analysis/all.preclustered.fasta"""
    cmd2 = """vsearch --uchime_deno v_analysis/all.preclustered.fasta \
	--sizein --sizeout --fasta_width 0 \
	--nonchimeras v_analysis/all.denovo.nonchimeras.fastaq"""
    cmd3 = """vsearch --uchime_ref v_analysis/all.denovo.nonchimeras.fastaq \
	--db /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/rdp_gold.fa --sizein --sizeout --fasta_width 0 \
	--nonchimeras v_analysis/all.ref.nonchimeras.fasta"""
    cmd4 = """perl /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/map.pl splited/derep.fa v_analysis/all.preclustered.uc v_analysis/all.ref.nonchimeras.fasta > v_analysis/all.nonchimeras.derep.fasta"""
    cmd5 = """perl /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/map.pl splited/derep.fa splited/all.derep.uc v_analysis/all.nonchimeras.derep.fasta > v_analysis/all.nonchimeras.fasta"""


def task_cluster():
    cmd = """vsearch --cluster_size v_analysis/all.nonchimeras.fasta --id 0.97 \
	--sizein --sizeout --fasta_width 0 \
	--uc v_analysis/all.clustered.uc \
	--relabel OTU --centroids v_analysis/all.otus.fasta"""
    cmd2 = """vsearch --usearch_global splited/filtered_uparsed.fa --db v_analysis/all.otus.fasta --strand plus --id 0.97 --uc v_analysis/map.txt --otutabout v_analysis/otu_raw.tab    """


def perform_vsearch_pipelines(vsearch_path,
                              ):
    vsearch_path = vsearch_path if vsearch_path else '/usr/bin/vsearch'
