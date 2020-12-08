import click
import classifycopynumber.classify
import classifycopynumber.parsers
import classifycopynumber.plots
import os

@click.command()
@click.argument('genes_gtf')
@click.argument('results_dir')
@click.argument('sample')
@click.argument('amps')
@click.argument('dels')

@click.option('--remixt_h5_filename', help='Remixt H5')
@click.option('--hmmcopy_csv_filename', help='hmmcopy csv')
@click.option('--plot', help='bool to generate plots')


def main(genes_gtf, results_dir, sample, amps, dels, remixt_h5_filename=None, hmmcopy_csv_filename=None, plot=False):
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    genes = classifycopynumber.parsers.read_gene_data(genes_gtf)

    if remixt_h5_filename is None and hmmcopy_csv_filename is None:
        raise click.ClickException('One of remixt_h5_filename, hmmcopy_csv_filename required')

    if remixt_h5_filename is not None:
        cn, stats = classifycopynumber.parsers.read_remixt(remixt_h5_filename)
        cn_cols = ['major_raw','minor_raw', 'copy','major_1',
        'minor_1','major_2','minor_2']
        ploidy = stats["ploidy"]

    elif hmmcopy_csv_filename is not None:
        cn, ploidy = classifycopynumber.parsers.read_hmmcopy(hmmcopy_csv_filename, sample)
        cn_cols=["state", "copy"]
    print(cn, cn.columns)
    genes_of_interest = classifycopynumber.parsers.compile_genes_of_interest(genes)

    gene_cn = classifycopynumber.classify.calculate_gene_cn_data(cn, genes_of_interest, cn_cols)

    amp_genes = classifycopynumber.classify.label_amplifications(gene_cn, ploidy, sample)
    hdel_genes = classifycopynumber.classify.label_deletions(gene_cn, ploidy, sample)

    amp_genes.to_csv(amps)
    hdel_genes.to_csv(dels)

    if plot == "True":
        plots_dir = os.path.join(results_dir, "plots")
        if not os.path.exists(plots_dir):
            os.makedirs(plots_dir)
        log_change_plot = os.path.join(plots_dir, "log_change_amplification_genes.png")
        classifycopynumber.plots.plot_log_change(amp_genes, log_change_plot)

        amp_genes = amp_genes[amp_genes.pass_filter==True]
        log_change_plot = os.path.join(plots_dir, "log_change_amplified_amplification_genes.png")
        classifycopynumber.plots.plot_log_change(amp_genes, log_change_plot)

        log_change_plot = os.path.join(plots_dir, "log_change_hdel_genes.png")
        classifycopynumber.plots.plot_log_change(hdel_genes, log_change_plot)

        hdel_genes = hdel_genes[hdel_genes.pass_filter==True]
        log_change_plot = os.path.join(plots_dir, "log_change_hdel_deleted_genes.png")
        classifycopynumber.plots.plot_log_change(hdel_genes, log_change_plot)
