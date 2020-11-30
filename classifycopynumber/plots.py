import seaborn
import matplotlib.pyplot as plt


def plot_log_change(amp_data, fig):
    if amp_data.empty:
        f, _ =plt.subplots()
        f.savefig(fig)
        return
    gene_order = amp_data.groupby(['gene_name'])['log_change'].mean().sort_values()
    g = seaborn.factorplot(x='log_change', y='gene_name', data=amp_data, order=gene_order.index, kind='box', aspect=0.75)
    g.axes[0][0].set_xlabel('Log2 change from ploidy')
    g.axes[0][0].set_ylabel('Gene')
    g.fig.savefig(fig, bbox_inches='tight')

#I think this is  really only useful for multple samples
# def plot_amp_matrix(amp_data, fig):
#     # amp matrix plot
#     plot_data =amp_data[amp_data.log_change>0.6].set_index(['gene_name'])['log_change']
#     # plot_data = amp_data.query('log_change > 1.').set_index(['gene_name'])['log_change']
#     plot_data.index.name = 'Gene'
#     print(plot_data)
#     # plot_data.columns.name = 'Patient'
#     # gene_order = plot_data.isnull().sum(axis=1).sort_values().index
#     # plot_data = plot_data.reindex(index=gene_order)
#     # sample_order = (plot_data.notnull() * 1).apply(lambda a: ''.join([str(b) for b in a])).sort_values(ascending=False)
#     # sample_order = list(sample_order.index)
#     # for s in cn_data.keys():
#     #     if s not in sample_order:
#     #         sample_order.append(s)
#     # plot_data = plot_data.reindex(columns=sample_order)

#     fig = plt.figure(figsize=(12, 6))
#     ax = seaborn.heatmap(plot_data, cmap='Reds', annot=True)
#     ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[1])))
#     ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[0])))
#     ax.xaxis.grid(True, which='minor', linestyle=':')
#     ax.yaxis.grid(True, which='minor', linestyle=':')
#     ax.spines['bottom'].set_visible(True)
#     ax.spines['bottom'].set_linestyle(':')
#     ax.spines['left'].set_visible(True)
#     ax.spines['top'].set_visible(True)
#     ax.spines['right'].set_visible(True)
#     fig.savefig(fig, bbox_inches='tight')


# def plot_hdels(hdel_data, fig):
#     # hdel plot
#     print(hdel_data)
#     plot_data = (
#         hdel_data
#         .groupby(['gene_name'])['hdel_width'].max()
#         .unstack())
#     plot_data = plot_data / 1000
#     plot_data.index.name = 'Gene'
#     plot_data.columns.name = 'Patient'
#     gene_order = plot_data.isnull().sum(axis=1).sort_values()
#     plot_data = plot_data.reindex(index=gene_order.index)
#     plot_data = plot_data.reindex(columns=sample_order)

#     fig = plt.figure(figsize=(12, 8))
#     ax = seaborn.heatmap(plot_data, cmap='Blues', annot=False)
#     ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[1])))
#     ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(plot_data.shape[0])))
#     ax.xaxis.grid(True, which='minor', linestyle=':')
#     ax.yaxis.grid(True, which='minor', linestyle=':')
#     ax.spines['bottom'].set_visible(True)
#     ax.spines['bottom'].set_linestyle(':')
#     ax.spines['left'].set_visible(True)
#     ax.spines['top'].set_visible(True)
#     ax.spines['right'].set_visible(True)
#     ax.tick_params(axis='y', rotation=0)
#     fig.savefig(fig, bbox_inches='tight')
