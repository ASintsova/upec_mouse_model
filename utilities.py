import configparser
import os
from scipy import stats


cols = ["#0084ff", "#44bec7", "#ffc300", "#fa3c4c", "#d696bb"]


def process_config(config_file="config"):

    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = {}
    for section in config.sections():
        config_dict[section] = {name: value for name, value in config.items(section)}
    return config_dict


def invnorm(x):
    return stats.norm.ppf((x.rank() - 0.5)/x.count())


def process_gff(gff_file, feat='CDS', id_sym="gene_id="):
    """
    Only been counting features that are 'CDS', consequently here also only looking at
    lines that have CDS in them
    :param gff_file:
    :return: dictionary of gene_id: gene length in kb
    """
    gene_to_gene_length = {}
    with open(gff_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                break
            elif line.startswith(("#", " ")) or len(line) == 0:
                continue
            elif line.split('\t')[2] != feat:
                continue
            else:
                start = int(line.split("\t")[3].strip())
                end = int(line.split("\t")[4].strip())
                gene_length = abs(end - start)/1000
                #prokka = line.split("\t")[-1].split(";")[0].strip(id_sym)
                prokka = line.split("\t")[-1].split(id_sym)[1].split(";")[0]
                # This would give me the prokka id
                gene_to_gene_length[prokka] = gene_length
    return gene_to_gene_length


def process_count_file(count_file):
    line = (l.split("\t") for l in open(count_file))
    counts = {g[0]: int(g[1].strip()) for g in line}
    return counts


def calculate_tpm(counts_dict, gene_len_dict):
    total_rpk = 0
    temp_rpk = {}
    for gene, count in counts_dict.items():
        if gene.startswith("__"):  # HTSeq specific: end of file has __total_mapped reads, etc.
            continue
        try:
            gene_length = gene_len_dict[gene]
        except KeyError:
            continue  # skipping genes we don't have length for
        else:
            rpk = count/gene_length
            total_rpk += rpk
            temp_rpk[gene] = rpk
    total_rpk /= 1000000  # Make sure this is a million
    tpm = {gene: rpk/total_rpk for gene, rpk in temp_rpk.items()}
    return tpm


def normalize_counts_to_tpm_one_file(cf, gff, feat='CDS', id_sym='gene_id='):
    counts_dict = process_count_file(cf)
    strain = os.path.basename(cf).split("_")[0]
    if os.path.isdir(gff):
        gff_file = os.path.join(gff, "{}.gff".format(strain))
    else:
        gff_file = gff
    gene_len_dict = process_gff(gff_file, feat, id_sym)
    tpm = calculate_tpm(counts_dict, gene_len_dict)
    return tpm


def normalize_counts_to_tpm(counts_dir, gff_dir, out_dir, feat='CDS', id_sym='gene_id='):
    """
    Assumes names of counts strats with strain, and gff named strain.gff
    :param counts_dir:
    :param gff_dir:
    :param out_dir:
    :return:
    """
    count_files = [os.path.join(counts_dir, f) for f in os.listdir(counts_dir)]
    all_tpms = {}
    for cf in count_files:
        print(cf)
        if "_counts" not in os.path.basename(cf):
            continue
        tpm = normalize_counts_to_tpm_one_file(cf, gff_dir, feat, id_sym)
        #out_file = "{}_tpm.csv".format(os.path.basename(cf))
        #out_path = os.path.join(out_dir, out_file)
        #with open(out_path, "w") as fo:
         #   for gene, t in tpm.items():
          #      fo.write("{},{}\n".format(gene, t))
        prefix = os.path.basename(cf).split("_trimmed")[0] # this would be specific to my naming convention
        all_tpms[prefix] = tpm
    return all_tpms



"""
Plotting PCA elipses:

__author__:

"""

def plot_point_cov(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


def find_pc1_pc2(df, meta):
    df = df.T
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(df)
    pDf = (pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
           .set_index(df.index))
    pc1_var = round(pca.explained_variance_ratio_[0] * 100, 2)
    pc2_var = round(pca.explained_variance_ratio_[1] * 100, 2)
    pDf2 = pDf.merge(meta, left_index=True, right_index=True)
    return pDf2, pc1_var, pc2_var


def plotPCA(pDf, pc1_var, pc2_var, colorby, col, nameby="", el=False):
    sns.set_style("ticks")
    sns.set_context("notebook", font_scale=2.2)
    group = pDf[colorby].unique()
    assert len(group) <= len(col)
    fig = plt.figure(figsize=(8, 8))
    for g, c in zip(group, col):
        df = pDf[pDf[colorby] == g]
        x, y = df[["PC1"]].values, df[["PC2"]].values
        ax = plt.scatter(x, y, c=c, s=150, label=g)
        if el:
            pts = np.asarray([[float(a), float(b)] for a, b in zip(x, y)])
            plot_point_cov(pts, nstd=2, alpha=0.1, color=c)
        if nameby:
            labels = df[nameby]
            for label, pc1, pc2 in zip(labels, x, y):
                plt.annotate(label, xy=(pc1, pc2), xytext=(-5, 7), textcoords="offset points",fontsize=14)
        plt.xlabel('Principal Component 1, {} %'.format(pc1_var), )
        plt.ylabel('Principal Component 2, {} %'.format(pc2_var), )
        #plt.xticks(fontsize=16)
        #plt.yticks(fontsize=16)
        plt.legend(frameon=True)
    return fig



