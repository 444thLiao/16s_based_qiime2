import re,gzip
from tqdm import tqdm
from default_params import *
from utils import get_files, data_parser
from Bio import SeqIO

def p_sto(col, stodge, is_id=False):
    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    for v in col:
        if not is_id:
            v = re.compile(''.join([iupac[_.upper()] for _ in v]))
        stodge.append(v)

def parse_metadata(metadata,
                   id_col='',
                   fb_col='',
                   rb_col='',
                   fp_col='',
                   rp_col='',
                   ):
    ids = []
    fb = []
    rb = []
    fp = []
    rp = []
    col2sto = dict(zip([id_col, fb_col, rb_col, fp_col, rp_col],
                       [ids, fb, rb, fp, rp]))
    metadf = data_parser(metadata, ft='csv')

    for col, sto in col2sto.items():
        if col and col in metadf.columns:
            column = metadf.loc[:, col]
            if col == id_col:
                p_sto(column, sto, is_id=True)
            else:
                p_sto(column, sto)
    return ids, fb, rb, fp, rp

def process_pair(read1_data,
                 read2_data,
                 forward_primers=None,
reverse_primers=None,
attempt_read_orientation=False,
output_fastq1=None,
output_fastq2=None,
                 ):

    # init
    is_reversed = False
    is_close_circle = False
    ## force re orientation
    for curr_primer in forward_primers:
        bc1_end = curr_primer.search(str(read1_data.seq))
        if bc1_end is not None:
            read1 = read1_data
            read2 = read2_data
            bc1_end = bc1_end.start
            break
        bc1_end = curr_primer.search(str(read2_data.seq))
        if bc1_end is not None:
            read1 = read2_data
            read2 = read1_data
            is_reversed = True
            bc1_end = bc1_end.start()
            break
    # Check reverse primers if forward primers not found
    for curr_primer in reverse_primers:
        bc2_end = curr_primer.search(str(read1_data.seq))
        if bc2_end is not None:
            read1 = read2_data
            read2 = read1_data
            bc2_end = bc2_end.start()
            if is_reversed:
                is_close_circle = True
            break
        bc2_end = curr_primer.search(str(read2_data.seq))
        if bc2_end is not None:
            read1 = read1_data
            read2 = read2_data
            bc2_end = bc2_end.start()
            if not is_reversed:
                is_close_circle = True
            break

    if is_close_circle and attempt_read_orientation:
        # 如果闭环,且要fix方向,则正常输出
        read1_output2file = output_fastq1
        read2_output2file = output_fastq2
    elif not attempt_read_orientation and is_close_circle:
        # 如果闭环,且不想fix方向,那么,按照原来的方向
        # 如果is_reversed,说明read1/2互换了,read1肯定到read1_output2file,所以此时,read1肯定到read1_output2file应该是
        #       output_fastq2
        read1_output2file = output_fastq2 if is_reversed else output_fastq1
        read2_output2file = output_fastq1 if is_reversed else output_fastq2
    else:
        # 如果不闭环,那么说明这个切都切不了..算了,告辞..不要了
        read1_output2file = None
        read2_output2file = None


    if bc1_end and bc2_end:  # self_add  by liaoth
        # print 'test successed'
        bc = read1[:bc1]
        bc_read1 = read1[bc1_end:]
        bc_read2 = read2[bc2_end - bc2_len:bc2_end]
    else:  # self_add  by liaoth
        bc_read1 = read1[0:bc1_len]
        bc_read2 = read2[0:bc2_len]


    if rev_comp_bc1:
        bc_read1 = str(DNA(bc_read1).rc())
        bc_qual1 = bc_qual1[::-1]
    if rev_comp_bc2:
        bc_read2 = str(DNA(bc_read2).rc())
        bc_qual2 = bc_qual2[::-1]

def process_single():
    pass

def demux_single(seqfile1,
               seqfile2,
               ):

    if seqfile2:
        if '.gz' in seqfile1:
            seqfile1 = gzip.open(seqfile1,'rt')
            seqfile2= gzip.open(seqfile2,'rt')

        for read1,read2 in zip(SeqIO.parse(seqfile1,format='fastq'),
                               SeqIO.parse(seqfile2,format='fastq')):
            process_pair(read1,
                         read2,
                         )
    else:
        process_single()