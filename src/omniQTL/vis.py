from .utils import *
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

class GeneTxPlot():
    def __init__(self):
        pass

    def read_bed12(self, bed12='Homo_sapiens.GRCh38.115.bed12'):
        self.df = pd.read_table(bed12, header=None, sep='\t', low_memory=False)
        self.df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                           'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                           'geneID', 'geneName', 'geneBiotype', 'transcriptID', 'transcriptName', 'transcriptBiotype', 'tag']
        return self.df

    def subset_gene(self, gene='GINS4', canonical_only=True, flank=1000, geneBiotype_include=['protein_coding'], geneName_exclude=['.']):
        wh = self.df['geneName'].isin([gene])
        df = self.df[wh]
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df['tag']])
            df = df[wh]

        chrom = df['chrom'].iloc[0]
        start = df['chromStart'].min()
        end = df['chromEnd'].max()
        window = [chrom, start - flank, end + flank]

        wh1 = (self.df['chrom'] == window[0]) & (self.df['chromStart'] > window[1]) & (self.df['chromEnd'] < window[2])
        wh2 = self.df['geneBiotype'].isin(geneBiotype_include)
        wh3 = ~self.df['geneName'].isin(geneName_exclude)
        df = self.df[wh1 & wh2 & wh3]
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df['tag']])
            df = df[wh]

        self.subset_file = gene + '.bed12'
        df.to_csv(self.subset_file, header=None, index=None, sep='\t')
        return window

    def plot_gene_tx(self, bed12_file='GINS4.bed12', window=['1', 116900000, 117000000], flank_no_window=10000, exon_height=0.2, label_shift=1000, fontsize=8, ax=None, show_x_ticks=True, cmap='pastel'):
        cmap = sns.color_palette(cmap)
        if ax is None:
            fig, ax = plt.subplots()

        df = self.read_bed12(bed12_file)

        L = df[['chromStart', 'chromEnd', 'transcriptName']].values
        G = self._find_non_overlapping_groups(L)

        if not window:
            chrom = str(df['chrom'].iloc[0])
            x_min = df['chromStart'].min() - flank_no_window
            x_max = df['chromEnd'].max() + flank_no_window
        else:
            chrom = str(window[0])
            x_min = window[1]
            x_max = window[2]
        if chrom.find('chr') == -1:
            chrom = 'chr' + chrom
        ax.set_xlim(x_min, x_max)

        if show_x_ticks:
            ax.set_xticks([x_min, int((x_min + x_max)/2), x_max])
            ax.set_xticklabels([x_min, chrom, x_max], fontsize=6)
        else:
            ax.set_xticks([])

        ax.set_yticks([])
        ax.set_ylim(0, len(set(G.values())) + 1)

        for n in range(df.shape[0]):
            start = df['chromStart'].values[n]
            end = df['chromEnd'].values[n]
            strand = df['strand'].values[n]
            strand = 0 if strand == '+' else 1
            tx_name = df['transcriptName'].values[n]
            gene_name = df['geneName'].values[n]
            blockCount = df['blockCount'].values[n]
            blockSizes = df['blockSizes'].values[n].split(',')
            blockStarts = df['blockStarts'].values[n].split(',')
            color = cmap[strand]
            y = G[tx_name]
            ax.plot([start, end], [y, y], color=color)
            for b in range(blockCount):
                block_start = start + int(blockStarts[b])
                block_size = int(blockSizes[b])
                ax.add_patch(Rectangle((block_start, y - exon_height/2), block_size, exon_height, color=color))
            ax.text(end + label_shift, y, gene_name, va='center', fontsize=fontsize)

    def _find_non_overlapping_groups(self, ranges=[(1, 3), (2, 5), (6, 8), (9, 10), (4, 7)]):
        # Sort the ranges by starting point (and by ending point in case of ties)
        ranges = sorted(ranges, key=lambda x: (x[0], x[1]))

        groups = []

        # Iterate through each range
        for range in ranges:
            placed = False
            # Try to place the current range in an existing group
            for group in groups:
                # Check if the current range overlaps with the last range in the current group
                if group[-1][1] < range[0]:
                    group.append(range)
                    placed = True
                    break
            # If the range couldn't be placed in any existing group, create a new group
            if not placed:
                groups.append([range])
        g = {}
        for i, group in enumerate(groups):
            for x in group:
                g[x[2]] = i + 1
        return(g)

class LocusZoomPlot(GeneTxPlot):
    def __init__(self):
        pass

    def scatter_plot(self, df, ax, window=['1', 116900000, 117000000], x='pos', y='pv', s=6, lw=0, color='C0', cmap='Set2', ylabel='-log10(p)', hue=None):
        ax.set_xlim([window[1], window[2]])
        if hue:
            sns.scatterplot(x=x, y=y, s=s, data=df, ax=ax, linewidth=lw, hue=hue, palette=cmap)
        else:
            sns.scatterplot(x=x, y=y, s=s, data=df, ax=ax, linewidth=lw, color=color)

        ax.set_xticks([])
        ax.set_xlabel('')
        ax.set_ylim(0, max(ax.get_ylim()[1], 10))
        ax.set_ylabel(ylabel)
        if ax.get_legend():
            ax.get_legend().remove()

    def bgzip_to_df(self, bgzip_file='', window=['1', 116900000, 117000000], gene='PTGFRN', extra_filter=None, adding_ld=True, var_id='', bfile=''}):
        header = pd.read_table(bgzip_file, header=0, nrows=0, sep='\t')
        tb = tabix.open(bgzip_file)
        chrom = window[0]
        if chrom.find('chr') == -1:
            chrom = 'chr' + chrom
        res = tb.query(chrom, window[1], window[2])
        df = pd.DataFrame(res)
        df.columns = header.columns
        if gene:
            wh = [True if gene in x.split('_')[-1].split(',') else False for x in df['phe_id']]
            df = df[wh]

        if gene and len(df['phe_id'].unique()) > 1:
            if extra_filter is None:
                # for caQTL, if there are multiple peaks of the same gene, just take the one with the smallest p-value if no extra_filter is provided
                L = []
                for gi, g in df.groupby('phe_id'):
                    L.append([g, g['nom_pval'].min()])
                L = sorted(L, key=lambda x: x[-1])
                df = L[0]
            else:
                # for caQTL, if there are multiple peaks of the same gene, take the one with the peak name containing extra_filter
                for gi, g in df.groupby('phe_id'):
                    if gi.find(extra_filter) != -1:
                        df = g
                        break

        if adding_ld:
            pass

if __name__ == '__main__':
    gene = 'PTGFRN'
    gtp = GeneTxPlot()
    gtp.read_bed12(bed12='Homo_sapiens.GRCh38.115.bed12')
    gtp.subset_gene(gene=gene)
    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot()
    gtp.plot_gene_tx(bed12_file=f'{gene}.bed12', ax=ax)
    lzp = LocusZoomPlot()

