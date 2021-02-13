import wgs_analysis.algorithms.cnv
import numpy as np
from classifycopynumber.transformations import calculate_log_change


def label_amplifications(data, ploidy):

    data=data[data.cn_type=="amplification"]
    data = data.merge(calculate_log_change(data, ploidy), on="gene_name")
    data["pass_filter"] =  data.log_change > 1
    return data


def label_deletions(data, ploidy):

    data=data[data.cn_type=="deletion"]
    data["pass_filter"] = data["copy"] < 0.5
    return data


def calculate_gene_cn_data(aggregated_cn_data, cgc_genes, copy_number_columns):
    return wgs_analysis.algorithms.cnv.calculate_gene_copy(aggregated_cn_data, cgc_genes,copy_number_columns)

    return genes_cn_data.reset_index()


    

