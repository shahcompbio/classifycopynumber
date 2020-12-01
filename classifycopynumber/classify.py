import wgs_analysis.algorithms.cnv
import numpy as np
from classifycopynumber.transformations import calculate_log_change


def label_amplifications(data, ploidy):

    data=data[data.cn_type=="amplification"]

    data = data.merge(calculate_log_change(data, ploidy), on="gene_name")
    # normalize = (
    #     data.groupby('gene_name')['overlap_width']
    #     .sum().rename('sum_overlap_width').reset_index())

    # data['total_raw_weighted'] = data['copy'] * data['overlap_width']

    # data = data.groupby(['gene_name'])['total_raw_weighted'].sum().reset_index()
    # data = data.merge(normalize)
    # data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
    # data['log_change'] = np.log2(data['total_raw_mean'])/ploidy
    data["pass_filter"] =  data.log_change > 1

    return data


def label_deletions(data, ploidy):

    data=data[data.cn_type=="deletion"]
    print(data["copy"])
    data = data.merge(calculate_log_change(data, ploidy), on="gene_name")
    print(data["copy"])
    data["pass_filter"] = data["copy"] < 0.5
    # data = data[data['copy'] < 0.5]
    
    # normalize = data.groupby(['gene_name'])['overlap_width'].sum().rename('sum_overlap_width').reset_index()

    # data['total_raw_weighted'] = data['copy'] * data['overlap_width']

    # data = data.groupby(['gene_name'])['total_raw_weighted'].sum().reset_index()
    # data = data.merge(normalize)
    # data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
    # data['log_change'] = np.log2(data['total_raw_mean'])/ploidy
    # data = data[data['sum_overlap_width'] > 10000]

    return data


# def label_hdels(genes_cn_data):

  
#         data = data[data['total_raw'] < 0.5]
#         data = data.groupby(['gene_id'])['overlap_width'].sum().rename('hdel_width').reset_index()
#         data = data[data['hdel_width'] > 10000]
#         data['sample'] = sample

#         hdel_data.append(data[['gene_id', 'hdel_width', 'sample']])

#     hdel_data = pd.concat(hdel_data)

#     gene_cols = [
#         'gene_id',
#         'chromosome',
#         'gene_start',
#         'gene_end',
#         'gene_name',
#         'cn_type',
#         'Tumour Types(Somatic)',
#     ]
#     hdel_data = hdel_data.merge(cgc_genes[gene_cols].query('cn_type == "deletion"').drop_duplicates())
#     return hdel_data
#     hdel_data.to_csv('../data/hdels.csv.gz', index=False)


def calculate_gene_cn_data(aggregated_cn_data, cgc_genes, copy_number_columns):
    return wgs_analysis.algorithms.cnv.calculate_gene_copy(aggregated_cn_data, cgc_genes,copy_number_columns)

    #not per cell
    # genes_cn_data= aggregated_cn_data.groupby(group_label_col).apply(lambda group: 
    #     wgs_analysis.algorithms.cnv.calculate_gene_copy(group, cgc_genes,copy_number_columns)
    # )
    return genes_cn_data.reset_index()


    

