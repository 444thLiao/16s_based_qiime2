import os
from os.path import join,dirname

import luigi

from luigi_pipelines import fileparser, run_cmd, config, valid_path, basic_luigi_task, visulize_seq, summarized_tasks,tabulate_seq
import sys
sys.path.insert(0, dirname(__file__))

class import_data(basic_luigi_task):
    def output(self):
        return luigi.LocalTarget(join(str(self.odir),
                                      config.preprocess_path,
                                      config.raw_seq_path + '.qza'), )

    def run(self):
        valid_path(self.output().path, check_ofile=True)
        df = fileparser(self.tab)
        manifest_path = df.generate_manifest(os.path.join(str(self.odir),
                                                          "data_manifest"))
        cmd = "qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path {manifest} \
  --output-path {prefix}.qza \
  --input-format PairedEndFastqManifestPhred33".format(
            manifest=manifest_path,
            prefix=self.output().path.replace('.qza', ''), )
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())


class joined_fastq(basic_luigi_task):

    def requires(self):
        return import_data(tab=self.tab,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path)

    def output(self):
        return luigi.LocalTarget(join(str(self.odir),
                                      config.preprocess_path,
                                      config.joined_seq_path + '.qza'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        extra_str = ''
        for p, val in config.join_params.items():
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)
        cmd = "qiime vsearch join-pairs --i-demultiplexed-seqs {input_file} --o-joined-sequences {ofile}".format(input_file=self.input().path,
                                                                                                                 ofile=self.output().path) + extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )


class qualityFilter(basic_luigi_task):
    def requires(self):
        return joined_fastq(tab=self.tab,
                            odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace(config.joined_seq_path,
                                                           config.joined_qc_seq_path))

    def run(self):
        extra_str = ''
        for p, val in config.qc_joined_params.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)
        cmd = "qiime quality-filter q-score-joined --i-demux {input_qza} --o-filtered-sequences {output_seq} --o-filter-stats {output_stats}".format(
            input_qza=self.input().path,
            output_seq=self.output().path,
            output_stats=self.output().path.replace(
                '.qza', '-stats.qza'))
        # path of stats used in `share_tasks.summaried_tasks`

        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )


class run_dada2(basic_luigi_task):
    mission = 'dada2'

    def requires(self):
        return import_data(tab=self.tab,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path)

    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                config.profiled_tab_path),
                           join(self.odir,
                                "%s_output" % self.mission,
                                config.representative_sequence_path),
                           join(self.odir,
                                "%s_output" % self.mission,
                                config.process_stats_path)]
                          ))
        return ofiles

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        extra_str = ''
        for p, val in config.dada2_args.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)

        cmd = """qiime dada2 denoise-paired --i-demultiplexed-seqs {input_file} --o-representative-sequences {rep_seq} --o-table {profiling_tab} --o-denoising-stats {stats_file}""".format(
            input_file=self.input().path,
            rep_seq=self.output()[1].path,
            profiling_tab=self.output()[0].path,
            stats_file=self.output()[2].path, )
        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class run_deblur(basic_luigi_task):
    mission = "deblur"

    def requires(self):
        return qualityFilter(tab=self.tab,
                             odir=self.odir,
                             dry_run=self.dry_run,
                             log_path=self.log_path)

    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                config.profiled_tab_path),
                           join(self.odir,
                                "%s_output" % self.mission,
                                config.representative_sequence_path),
                           join(self.odir,
                                "%s_output" % self.mission,
                                config.process_stats_path)]
                          ))
        return ofiles

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        extra_str = ''
        for p, val in config.deblur_args.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)

        cmd = """qiime deblur denoise-16S --i-demultiplexed-seqs {input_file} --o-representative-sequences {rep_seq} --o-table {profiling_tab} --o-stats {stats_file}""".format(
            input_file=self.input().path,
            rep_seq=self.output()[1].path,
            profiling_tab=self.output()[0].path,
            stats_file=self.output()[2].path, )
        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


############################################################
class joined_summarize(visulize_seq):
    def requires(self):
        return joined_fastq(tab=self.tab,
                            odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path)


class view_rep_seq(tabulate_seq):
    analysis_mission = luigi.Parameter(default="both")
    def requires(self):
        kwargs = dict(tab=self.tab,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        if self.analysis_mission == "dada2":
            return run_dada2(**kwargs)
        elif self.analysis_mission == "deblur":
            return run_deblur(**kwargs)
        else:
            raise Exception("wrong analysis_mission %s" % self.analysis_mission)


# class dada2_summarize(visulize_table):
#     def requires(self):
#         pass
#     def output(self):
#         pass
#     def run(self):
#         pass


class summarize_all(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    log_path = luigi.Parameter(default=None)
    analysis_mission = luigi.Parameter(default="both")

    def requires(self):
        kwargs = dict(tab=self.tab,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        am = str(self.analysis_mission)
        required_tasks = {}
        if am == "both":
            required_tasks["deblur"] = run_deblur(**kwargs)
            required_tasks["dada2"] = run_dada2(**kwargs)
            required_tasks["view_rep_seq1"] = view_rep_seq(analysis_mission="dada2",
                                                          **kwargs)
            required_tasks["view_rep_seq2"] = view_rep_seq(analysis_mission="deblur",
                                                          **kwargs)
        elif am == "dada2":
            required_tasks["dada2"] = run_dada2(**kwargs)
            required_tasks["view_rep_seq"] = view_rep_seq(analysis_mission=am,
                                                          **kwargs)
        elif am == "deblur":
            required_tasks["deblur"] = run_deblur(**kwargs)
            required_tasks["view_rep_seq"] = view_rep_seq(analysis_mission=am,
                                                          **kwargs)
        else:
            raise Exception("wrong analysis_mission %s" % self.analysis_mission)
        required_tasks["joined"] = joined_summarize(**kwargs)
        required_tasks["joined_qc"] = qualityFilter(**kwargs)
        return required_tasks

    def output(self):
        ofile = os.path.join(self.odir,config.summarized_stats_path)
        return luigi.LocalTarget(ofile)

    def run(self):
        merged_df = {}
        for k, v in self.input().items():
            df = summarized_tasks(k,v)
            if df is not None:
                merged_df[k] = df
        import pandas as pd
        final_df = pd.concat([_ for _ in merged_df.values()],axis=1)
        final_df.to_csv(self.output().path)


############################################################
class get_tree(basic_luigi_task):
    """
    todo: do it myself...
    """
    prefix = luigi.Parameter()

    def requires(self):
        pass
        # return run_deblur(tab=self.tab,
        #                      odir=self.odir,
        #                      dry_run=self.dry_run,
        #                      log_path=self.log_path)
        # or dada2

    def output(self):
        return luigi.LocalTarget(config.root_tree_path)

    def run(self):
        # cmd = "qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {rep_seq} --o-alignment {aligned_seq} --o-masked-alignment {masked_aligned_seq} --o-tree {out_tree} --o-rooted-tree {out_rooted_tree}".format(
        #     rep_seq=self.input().path,
        #     aligned_seq='',
        #     masked_aligned_seq='',
        #     out_tree='',
        #     out_rooted_tree='',
        # )
        pass


if __name__ == '__main__':
    luigi.run()
