"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
import os
import pandas as pd
from utils import valid_path


class fileparser():
    def __init__(self, filename):
        filename = os.path.abspath(filename)
        self.df = pd.read_csv(filename, index_col=None,dtype=str)
        self.cols, self.df = validate_df(self.df, filename)
        self.df = self.df.set_index("sample_ID")

    def get_attr(self, col):
        if col == self.df.index.name:
            return list(self.df.index)
        if col not in self.cols:
            raise Exception("attr %s not in input df" % col)
        else:
            return list(self.df[col])

    def generate_manifest(self, opath):
        valid_path(opath, check_ofile=1)
        template_text = "sample-id,absolute-filepath,direction\n"
        r1_files = self.get_attr("path_R1")
        r2_files = self.get_attr("path_R2")
        ids = self.get_attr("sample_ID")
        for r1, r2, sid in zip(r1_files, r2_files, ids):
            template_text += ','.join([sid, r1, 'forward']) + '\n'
            template_text += ','.join([sid, r2, 'reverse']) + '\n'
        with open(opath, 'w') as f1:
            f1.write(template_text)
        return opath


def validate_df(df, filename):
    template_file = os.path.join(os.path.dirname(__file__),
                                 "data_input.template")
    columns_values = open(template_file).read().strip('\n').split(',')

    if set(df.columns) != set(columns_values):
        raise Exception("INPUT file has unknown header")

    if df["sample_ID"].duplicated().any():
        raise Exception("sample_ID has duplicated.")

    chdir = os.path.dirname(os.path.abspath(filename))
    for idx, row in df.iterrows():
        # auto implement filepath
        # so easy~~~
        df.loc[idx, "path_R1"] = os.path.join(chdir, row["path_R1"])
        df.loc[idx, "path_R2"] = os.path.join(chdir, row["path_R2"])
    return columns_values, df
