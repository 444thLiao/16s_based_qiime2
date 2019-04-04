from qiime2.plugins.deblur.methods import denoise_16S
from qiime2.plugins.vsearch.methods import uchime_denovo
from qiime2.plugins.vsearch.methods import cluster_features_de_novo


def OTU_pipelines(deblur_input,
table=table,
                             perc_identity=perc_identity,
                             threads=threads
                  ):
    cluster_features_de_novo(sequences=deblur_input,
                             table=table,
                             perc_identity=perc_identity,
                             threads=threads,)