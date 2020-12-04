import wgs_analysis.algorithms.cnv
import numpy as np
from classifycopynumber.transformations import calculate_log_change


def label_amplifications(data, ploidy, sample):

    data=data[data.cn_type=="amplification"]
    data.to_csv("/juno/work/shah/abramsd/RESULTS/dlp_cohort_oncoplot/DATA/copynumber_classification/signatures/{}_aplications_pre_log_calc.tsv".format(sample))
    data = data.merge(calculate_log_change(data, ploidy), on="gene_name")
    data["pass_filter"] =  data.log_change > 1
    return data


def label_deletions(data, ploidy, sample):

    data=data[data.cn_type=="deletion"]
    data.to_csv("/juno/work/shah/abramsd/RESULTS/dlp_cohort_oncoplot/DATA/copynumber_classification/signatures/{}_deletions_pre_log_calc.tsv".format(sample))
    data["pass_filter"] = data["copy"] < 0.5
    return data


def calculate_gene_cn_data(aggregated_cn_data, cgc_genes, copy_number_columns):
    return wgs_analysis.algorithms.cnv.calculate_gene_copy(aggregated_cn_data, cgc_genes,copy_number_columns)

    return genes_cn_data.reset_index()


    

