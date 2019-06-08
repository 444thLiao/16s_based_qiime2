import os

import luigi

from luigi_pipelines import run_cmd, config
from qiime2 import Artifact
import pandas as pd
import qiime2

class basic_luigi_task(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    log_path = luigi.Parameter(default=None)

    def get_log_path(self):
        if self.log_path is not None:
            return self.log_path
        else:
            return os.path.join(str(self.odir),
                                      str(self.log_path))


class visulize_seq(basic_luigi_task):
    """
    mainly for visualizing SingleFastqFormat qza
    """
    def output(self):
        return luigi.LocalTarget(self.input().path.replace(".qza",".qzv"))

    def run(self):
        cmd = "qiime demux summarize --i-data {input_f} --o-visualization {output_f}".format(
            input_f=self.input().path,
            output_f=self.output().path)
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )


class tabulate_seq(basic_luigi_task):
    """
    mainly for visualizing representative sequence
    """
    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace(".qza",".qzv"))

    def run(self):
        cmd = "qiime feature-table tabulate-seqs --i-data {input_f} --o-visualization {output_f}".format(
            input_f=self.input()[1].path,
            output_f=self.output().path)
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )

#
# class visulize_table(basic_luigi_task):
#     """
#     mainly for visualizing SingleFastqFormat qza
#     """
#     def output(self):
#         return luigi.LocalTarget(self.input().path.replace(".qza",".qzv"))
#
#     def run(self):
#         cmd = "qiime metadata tabulate --m-input-file {input_f} --o-visualization {output_f}".format(
#             input_f=self.input().path,
#             output_f=self.output().path)
#         run_cmd(cmd,
#                 dry_run=self.dry_run,
#                 log_file=self.get_log_path())
#         if self.dry_run:
#             run_cmd("touch %s" % self.output().path, dry_run=False, )


def summarized_tasks(key,values):
    if key == "dada2":
        stats_p = [_.path for _ in values if config.process_stats_path in _.path]
        if stats_p:
            artifact = Artifact.load(stats_p[0])
            dada2_df = artifact.view(qiime2.Metadata).to_dataframe()
            dada2_df.index = dada2_df.index.astype(str)
            dada2_df.append(pd.DataFrame([[]]))
            return dada2_df

    elif key == "deblur":
        stats_p = [_.path for _ in values if config.process_stats_path in _.path]
        if stats_p:
            artifact = Artifact.load(stats_p[0])
            deblur_df = artifact.view(pd.DataFrame)
            deblur_df.index = deblur_df.index.astype(str)
            return deblur_df
    elif key == "joined_qc":
        file_p = values.path.replace('.qza', '-stats.qza')
        # qualityFilter output stats path
        artifact = Artifact.load(file_p)
        joined_qc_df = artifact.view(pd.DataFrame)
        joined_qc_df.index = joined_qc_df.index.astype(str)
        return joined_qc_df