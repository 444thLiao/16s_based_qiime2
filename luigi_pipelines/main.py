import sys
from os.path import dirname
import luigi

sys.path.insert(0, dirname(dirname(__file__)))
sys.path.insert(0, dirname(__file__))



class main(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="OTU")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)
    def requires(self):
        pass
    def run(self):
        pass
