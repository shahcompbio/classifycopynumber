import click
import classifycopynumber.classify
import classifycopynumber.parsers
import classifycopynumber.plots
import wgs_analysis.algorithms.cnv
import os
import pandas as pd
import yaml


@click.command()
@click.argument('genes_gtf')
@click.argument('results_dir')
@click.argument('amps')
@click.argument('dels')
@click.option('--remixt_h5_filename', help='Remixt H5')
@click.option('--remixt_parsed_csv', help='Remixt H5')
@click.option('--hmmcopy_csv_filenames', help='hmmcopy csv',  multiple=True)
@click.option('--sample_ids', help='sample ids',  multiple=True)
@click.option('--plot', is_flag=True, help='bool to generate plots')
def main(genes_gtf, results_dir, amps, dels, remixt_h5_filename=None, remixt_parsed_csv=None, hmmcopy_csv_filenames=(), sample_ids=(), plot=False):
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    genes = classifycopynumber.parsers.read_gene_data(genes_gtf)
    genes_of_interest = classifycopynumber.parsers.compile_genes_of_interest(genes)

    if remixt_h5_filename is None and len(hmmcopy_csv_filenames) == 0 and remixt_parsed_csv is None:
        raise click.ClickException('One of remixt_h5_filename, hmmcopy_csv_filename required')

    if remixt_parsed_csv is not None:
        cn, stats = classifycopynumber.parsers.read_remixt_parsed_csv(remixt_parsed_csv)
        ploidy = stats['ploidy']

    if remixt_h5_filename is not None:
        cn, stats = classifycopynumber.parsers.read_remixt(remixt_h5_filename)
        ploidy = stats['ploidy']

    elif len(hmmcopy_csv_filenames) > 0:
        if len(sample_ids) == 0:
            raise click.ClickException('sample_ids required if providing hmmcopy')
        cn, ploidy = classifycopynumber.parsers.read_hmmcopy_files(hmmcopy_csv_filenames, sample_ids=sample_ids)

    gene_cn = wgs_analysis.algorithms.cnv.calculate_gene_copy(cn, genes_of_interest, ['total_raw'])

    amp_genes = classifycopynumber.classify.label_amplifications(gene_cn, ploidy)
    hdel_genes = classifycopynumber.classify.label_deletions(gene_cn, ploidy)

    amp_genes.to_csv(amps, index=False)
    hdel_genes.to_csv(dels, index=False)

    if plot:
        plots_dir = os.path.join(results_dir, 'plots')
        if not os.path.exists(plots_dir):
            os.makedirs(plots_dir)
        log_change_plot = os.path.join(plots_dir, 'log_change_amplification_genes.png')
        classifycopynumber.plots.plot_log_change(amp_genes, log_change_plot)

        amp_genes = amp_genes[amp_genes.pass_filter==True]
        log_change_plot = os.path.join(plots_dir, 'log_change_amplified_amplification_genes.png')
        classifycopynumber.plots.plot_log_change(amp_genes, log_change_plot)

        log_change_plot = os.path.join(plots_dir, 'log_change_hdel_genes.png')
        classifycopynumber.plots.plot_log_change(hdel_genes, log_change_plot)

        hdel_genes = hdel_genes[hdel_genes.pass_filter==True]
        log_change_plot = os.path.join(plots_dir, 'log_change_hdel_deleted_genes.png')
        classifycopynumber.plots.plot_log_change(hdel_genes, log_change_plot)
