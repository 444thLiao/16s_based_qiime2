from pp.pipelines import selective_p, tax_assign_qiime2, g_tree
from utils import *
from default_params import *

def preprocess():
    write_manifest(indir,
                   opath,
                   r1_format,
                   idpattern)
    # 准备序列的输入
    raw_seq = import_data_with_manifest(opath)
    # 将序列信息导入qiime的环境,可另存为qza

    raw_seq_eval_vis = seq_eval(raw_seq,
                                n=n)
    raw_seq_eval_vis.save(os.path.join(odir,
                                       raw_seq_vis_path))

    joined_seq, joined_seq_eval_vis, \
    joined_qc_seq, joined_qc_eval_vis, joined_qc_stats = join_seqs(raw_seq,
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
    raw_seq.save(os.path.join(odir, 'raw_data'))
    joined_qc_seq.save(os.path.join(odir,
                                    joined_seq_vis_path))
    joined_qc_eval_vis.save(os.path.join(odir,
                                         joined_qc_seq_vis_path))
    joined_seq.save(os.path.join(odir,
                                 joined_qc_seq_vis_path))
    joined_seq_eval_vis.save(os.path.join(odir,
                                          joined_seq_vis_path))

    joined_qc_stats.view(pd.DataFrame).to_csv(os.path.join(odir,
                                                           joined_qc_stats_tab_path), index=True)

    return raw_seq, joined_qc_seq


def run_pipelines(p, pipelines_args):
    pre_ = 'sOTU'
    tab, rep, tab_vis, seq_vis, stats_df = selective_p(p, pipelines_args)

    p_tab = tab.view(pd.DataFrame)
    name_dict = dict(zip(p_tab.columns, [pre_ + str(_) for _ in range(p_tab.shape[1])]))
    p_tab.columns = [name_dict[_] for _ in p_tab.columns]

    seq_f = glob(str(rep._archiver.data_dir) + '/*')
    if len(seq_f) != 1:
        raise Exception
    seq_f = seq_f[0]

    # output part
    tab.save(os.path.join(odir,
                              profiled_tab_path.format(prefix=p)))
    p_tab.to_csv(os.path.join(odir,
                              profiled_tab_path.format(prefix=p)),
                 sep='\t' if profiled_tab_path.endswith('.tab') else ',',
                 index=True)
    mv_seq(seq_f,
           opath=os.path.join(odir,
                              representative_sequence_path.format(prefix=p)),
           name_dict=name_dict)
    rep.save(os.path.join(odir,
                              representative_sequence_path.format(prefix=p)))
    tab_vis.save(os.path.join(odir,
                              profiled_tab_vis_path.format(prefix=p)))
    seq_vis.save(os.path.join(odir,
                              representative_sequence_vis_path.format(prefix=p)))
    stats_df.to_csv(os.path.join(odir,
                                 process_stats_path.format(prefix=p)))

    return tab, p_tab, rep


def after_otu(args):
    # assigning taxonomy, perform alpha,beta diversity
    tax_tab = tax_assign_qiime2(**args)
    rooted_tree = g_tree(**args)


if __name__ == '__main__':
    from utils import parse_param
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pipelines", help="Which kinds of pipelines you want to perform. \
                             [deblur|dada2]",
                        type=str, choices=['deblur', 'dada2'])
    parser.add_argument("--parameter", help="input file contains all parameters, template is place at the %s. This is a python script actually, so just use python code to write your code." % os.path.join(os.path.dirname(os.path.abspath(__file__)),'param.template'),default=os.path.join(os.path.dirname(__file__),'param.template'))

    args = parser.parse_args()
    p = args.pipelines
    parameter_f = args.parameter

    os.makedirs(odir,exist_ok=True)
    parse_param(parameter_f)
    ## 跑命令
    raw_seq, joined_qc_seq = preprocess()
    print("预处理完成,完成原始序列评估 与 joined, 去污染,去chimera,fix orientation")
    pipelines_args['deblur_input'] = joined_qc_seq
    pipelines_args['dada2_input'] = raw_seq
    tab, p_tab, rep = run_pipelines(p, pipelines_args)
    # after_otu_args['rep'] = rep
    # after_otu(after_otu_args)
