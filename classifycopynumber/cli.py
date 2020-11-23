import click
import classifycopynumber.classify
import classifycopynumber.parsers


@click.command()
@click.argument('genes_gtf')
@click.option('--remixt_h5_filename', help='Remixt H5')
@click.option('--hmmcopy_csv_filename', help='Remixt H5')
def main(genes_gtf, remixt_h5_filename=None, hmmcopy_csv_filename=None):
    genes = classifycopynumber.parsers.read_gene_data(genes_gtf)

    if remixt_h5_filename is None and hmmcopy_csv_filename is None:
        raise click.ClickException('One of remixt_h5_filename, hmmcopy_csv_filename required')
    
    if remixt_h5_filename is not None:
        cn, _ = classifycopynumber.parsers.read_remixt(remixt_h5_filename)

    elif hmmcopy_csv_filename is not None:
        cn, _ = classifycopynumber.parsers.read_hmmcopy(hmmcopy_csv_filename)

    from IPython import embed; embed(); raise


## Below has to be integrated....
# metadata paths are wrong: use package data instead
# change output paths

    amp_genes = pd.read_csv('../../metadata/Census_ampThu Apr 16 15_35_36 2020.csv')
    amp_genes['cn_type'] = 'amplification'

    del_genes = pd.read_csv('../../metadata/Census_delsThu Apr 16 15_36_24 2020.csv')
    del_genes['cn_type'] = 'deletion'

    cgc_genes = pd.concat([
        amp_genes[['Gene Symbol', 'cn_type', 'Tumour Types(Somatic)']],
        del_genes[['Gene Symbol', 'cn_type', 'Tumour Types(Somatic)']]], ignore_index=False)
    cgc_genes = cgc_genes.rename(columns={'Gene Symbol': 'gene_name'})

    gene_rename = {
        'NSD3': 'WHSC1L1',
        'MRE11': 'MRE11A',
        'SEM1': 'SHFM1',
    }

    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'MECOM', 'cn_type': 'amplification'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'PIK3CA', 'cn_type': 'amplification'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'RAD52', 'cn_type': 'amplification'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'RAD51', 'cn_type': 'deletion'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'KRAS', 'cn_type': 'amplification'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'RAD21', 'cn_type': 'amplification'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'CDK12', 'cn_type': 'amplification'}), ignore_index=True)

    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'SLFN11', 'cn_type': 'deletion'}), ignore_index=True)

    # From https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_HOMOLOGOUS_RECOMBINATION
    hr_genes = pd.read_csv('../../metadata/hr_genes.txt', sep='\t')

    for gene_name in hr_genes['Gene_Name']:
        cgc_genes = cgc_genes.append(pd.Series({'gene_name': gene_name, 'cn_type': 'deletion'}), ignore_index=True)

    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'STK19', 'cn_type': 'deletion'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'ELOF1', 'cn_type': 'deletion'}), ignore_index=True)
    cgc_genes = cgc_genes.append(pd.Series({'gene_name': 'ERCC6L2', 'cn_type': 'deletion'}), ignore_index=True)

    antigen_present_genes = [
        'CIITA',
        'IRF1',
        'PSME1',
        'PSME2',
        'PSME3',
        'ERAP1',
        'ERAP2',
        'PSMA7',
        'HSPBP1',
        'TAP1',
        'TAP2',
        'TAPBP',
        'CALR',
        'CANX',
        'PDIA3',
        'B2M',
    ]

    for gene_name in antigen_present_genes:
        cgc_genes = cgc_genes.append(pd.Series({'gene_name': gene_name, 'cn_type': 'deletion'}), ignore_index=True)

    cgc_genes['gene_name'] = cgc_genes['gene_name'].apply(lambda a: gene_rename.get(a, a))

    cgc_genes = cgc_genes.merge(genes, how='left')

    assert not cgc_genes['gene_start'].isnull().any()


    ##
    # Aggegate cn data
    ##

    aggregated_cn_data = {}

    for sample, cn in cn_data.items():
        aggregated_cn_data[sample] = wgs_analysis.algorithms.cnv.aggregate_adjacent(
            cn,
            value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            length_normalized_cols=['major_raw', 'minor_raw'],
        )


    ##
    # Calculate gene cn data
    ##

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


    ##
    # Calculate amplifications
    ##

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

    amp_data.to_csv('../data/amps.csv.gz', index=False)


    ##
    # Calculate hdels
    ##

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

    hdel_data.to_csv('../data/hdels.csv.gz', index=False)


    ##
    # Plots
    ##

    # gene log change plot

    gene_order = amp_data.groupby(['gene_name'])['log_change'].mean().sort_values()
    g = seaborn.factorplot(x='log_change', y='gene_name', data=amp_data, order=gene_order.index, kind='box', aspect=0.75)
    g.axes[0][0].set_xlabel('Log2 change from ploidy')
    g.axes[0][0].set_ylabel('Gene')
    g.fig.savefig('../plots/logchange.pdf', bbox_inches='tight')

    # amp matrix plot

    plot_data = amp_data.query('log_change > 1.').set_index(['gene_name', 'sample'])['log_change'].unstack()
    plot_data.index.name = 'Gene'
    plot_data.columns.name = 'Patient'
    gene_order = plot_data.isnull().sum(axis=1).sort_values().index
    plot_data = plot_data.reindex(index=gene_order)
    sample_order = (plot_data.notnull() * 1).apply(lambda a: ''.join([str(b) for b in a])).sort_values(ascending=False)
    sample_order = list(sample_order.index)
    for s in cn_data.keys():
        if s not in sample_order:
            sample_order.append(s)
    plot_data = plot_data.reindex(columns=sample_order)

    fig = plt.figure(figsize=(12, 6))
    ax = seaborn.heatmap(plot_data, cmap='Reds', annot=True)
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[1])))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[0])))
    ax.xaxis.grid(True, which='minor', linestyle=':')
    ax.yaxis.grid(True, which='minor', linestyle=':')
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linestyle(':')
    ax.spines['left'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    fig.savefig('../plots/amps.pdf', bbox_inches='tight')

    # hdel plot

    plot_data = (
        hdel_data
        .groupby(['gene_name', 'sample'])['hdel_width'].max()
        .unstack())
    plot_data = plot_data / 1000
    plot_data.index.name = 'Gene'
    plot_data.columns.name = 'Patient'
    gene_order = plot_data.isnull().sum(axis=1).sort_values()
    plot_data = plot_data.reindex(index=gene_order.index)
    plot_data = plot_data.reindex(columns=sample_order)

    fig = plt.figure(figsize=(12, 8))
    ax = seaborn.heatmap(plot_data, cmap='Blues', annot=False)
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[1])))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[0])))
    ax.xaxis.grid(True, which='minor', linestyle=':')
    ax.yaxis.grid(True, which='minor', linestyle=':')
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linestyle(':')
    ax.spines['left'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.tick_params(axis='y', rotation=0)
    fig.savefig('../plots/hdels.pdf', bbox_inches='tight')
