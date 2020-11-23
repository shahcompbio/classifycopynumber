
def label_amplifications(genes_cn_data):
    
    amp_data = []

    for sample, data in genes_cn_data.items():
        normalize = (
            data
            .groupby('gene_id')['overlap_width']
            .sum().rename('sum_overlap_width').reset_index())

        data['total_raw_weighted'] = data['total_raw'] * data['overlap_width']

        data = data.groupby(['gene_id'])['total_raw_weighted'].sum().reset_index()
        data = data.merge(normalize)
        data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
        data['sample'] = sample

        amp_data.append(data[['gene_id', 'total_raw_mean', 'sample']])

    amp_data = pd.concat(amp_data)

    gene_cols = [
        'gene_id',
        'chromosome',
        'gene_start',
        'gene_end',
        'gene_name',
        'cn_type',
        'Tumour Types(Somatic)',
    ]
    amp_data = amp_data.merge(cgc_genes[gene_cols].query('cn_type == "amplification"').drop_duplicates())
    amp_data = amp_data.merge(stats_data[['sample', 'ploidy']])
    amp_data['log_change'] = np.log2(amp_data['total_raw_mean'] / amp_data['ploidy'])

    return amp_data
    # amp_data.to_csv('../data/amps.csv.gz', index=False)


def label_hdels(genes_cn_data):

    hdel_data = []

    for sample, data in genes_cn_data.items():
        data = data[data['total_raw'] < 0.5]
        data = data.groupby(['gene_id'])['overlap_width'].sum().rename('hdel_width').reset_index()
        data = data[data['hdel_width'] > 10000]
        data['sample'] = sample

        hdel_data.append(data[['gene_id', 'hdel_width', 'sample']])

    hdel_data = pd.concat(hdel_data)

    gene_cols = [
        'gene_id',
        'chromosome',
        'gene_start',
        'gene_end',
        'gene_name',
        'cn_type',
        'Tumour Types(Somatic)',
    ]
    hdel_data = hdel_data.merge(cgc_genes[gene_cols].query('cn_type == "deletion"').drop_duplicates())
    return hdel_data
    hdel_data.to_csv('../data/hdels.csv.gz', index=False)


def calculate_gene_cn_data(aggregated_cn_data, cgc_genes):
    genes_cn_data = {}

    for sample, cn in aggregated_cn_data.items():
        cn['width'] = cn['end'] - cn['start']
        cn['total_raw'] = cn['major_raw'] + cn['minor_raw']

        genes_cn_data[sample] = wgs_analysis.algorithms.cnv.calculate_gene_copy(
            cn, cgc_genes,
            [
                'major_raw',
                'minor_raw',
                'total_raw',
                'major_1',
                'minor_1',
                'major_2',
                'minor_2',
            ])
    return genes_cn_data

    
    
def aggregate_cn_data(cn_data):
    aggregated_cn_data = {}

    for sample, cn in cn_data.items():
        aggregated_cn_data[sample] = wgs_analysis.algorithms.cnv.aggregate_adjacent(
            cn,
            value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            length_normalized_cols=['major_raw', 'minor_raw'],
        )
    return aggregated_cn_data
