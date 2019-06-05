import os
from os.path import join

import luigi

from luigi_pipelines import run_cmd, config

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


class visulize_table(basic_luigi_task):
    """
    mainly for visualizing SingleFastqFormat qza
    """
    def output(self):
        return luigi.LocalTarget(self.input().path.replace(".qza",".qzv"))

    def run(self):
        cmd = "qiime metadata tabulate --m-input-file {input_f} --o-visualization {output_f}".format(
            input_f=self.input().path,
            output_f=self.output().path)
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )