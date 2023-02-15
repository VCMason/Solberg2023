
def plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, prefix, color1='blue', label1='F - R'):
    import matplotlib.pyplot as plt
    import os

    path, plottitle = os.path.split(prefix)

    # dnorm is a dictionary
    # dnorm key is a float representig midpoint of range ex: 0.025 for 0.00-0.05
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse

    x = halfbins
    y1 = [dnorm[k][0] - dnorm[k][1] for k in halfbins]  # dnorm[k][1] # make reverse reads negative
    print(x)
    print(y1)
    plt.plot(x, y1, color=color1, label=label1)  # 'F'- 'R' or 'S' - 'A'
    plt.grid()
    plt.legend(loc='upper left')
    plt.title(plottitle)

    outpath = prefix + '.perbase.strandedness.AllClusters.Normalized.plots.pdf'
    print(outpath)
    plt.savefig(outpath)
    plt.show()
    plt.close()

    return


def plot_one_normalized_per_base_strandedness(dnorm, halfbins, prefix, color1='blue', color2='orange', label1='F', label2='R'):
    import matplotlib.pyplot as plt
    import os

    path, plottitle = os.path.split(prefix)

    # dnorm is a dictionary
    # dnorm key is a float representig midpoint of range ex: 0.025 for 0.00-0.05
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse

    x = halfbins
    y1 = [dnorm[k][0] for k in halfbins]
    y2 = [dnorm[k][1] for k in halfbins]
    y2 = [-1*number for number in y2]  # make negative
    print(x)
    print(y1)
    print(y2)
    plt.plot(x, y1, color=color1, label=label1)  # 'F' or 'S'
    plt.plot(x, y2, color=color2, label=label2)  # # 'R' or 'A'
    plt.grid()
    plt.legend(loc='upper left')
    plt.title(plottitle)

    outpath = prefix + '.perbase.strandedness.AllClusters.Normalized.plots.pdf'
    print(outpath)
    plt.savefig(outpath)
    plt.show()
    plt.close()

    return


def make_sns_pair_plot(matrix, columnlabels, pcaindex, prefix):
    import pandas as pd
    import matplotlib.pyplot as plt
    # Seaborn visualization library
    import seaborn as sns
    from sklearn.preprocessing import StandardScaler

    m = []
    for column in matrix:
        col = []
        for val in column:
            col.append(round(float(val), 2))
        m.append(col)

    dfx = pd.DataFrame(data=m, columns=columnlabels, index=pcaindex)
    pairdf = dfx.assign(mRNA=pcaindex)
    print(pairdf.head())
    # Create the default pairplot
    fig = sns.pairplot(pairdf, hue='mRNA')  # , hue='mRNA' # vars=pairdf.columns[:-1],
    outpath = prefix + '.rawvalues.pairplot.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    ### plot standardized values
    scaler = StandardScaler()
    scaler.fit(matrix)
    scaledmatrix = scaler.transform(matrix)
    dfx = pd.DataFrame(data=scaledmatrix, columns=columnlabels, index=pcaindex)
    #pcaindex = [int(i) for i in pcaindex]
    # pcaindex is list of '0' and '1's to specify hue
    pairdf = dfx.assign(mRNA=pcaindex)
    print(pairdf.head())
    # Create the default pairplot
    fig = sns.pairplot(pairdf, hue='mRNA')  # , hue='mRNA' # vars=pairdf.columns[:-1],
    outpath = prefix + '.pairplot.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    return


def make_pca(matrix, columnlabels, pcaindex, prefix, varexp=0.90):
    # matrix is list of lists, first column is rowIDs, first row is column names
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn import datasets
    from pca import pca
    from sklearn.preprocessing import StandardScaler

    ## iris = datasets.load_iris()
    ## label = iris.feature_names
    ## y = iris.target
    ## x = iris.data
    ## x = matrix
    #scaler = StandardScaler()
    #scaler.fit(matrix)
    #scaledmatrix = scaler.transform(matrix)
    ## matrix is list of lists
    ## columnlabels are the column names for the matrix
    ## pcaindex will be color of dots, i.e. if cluster is covered by >= number of mRNA
    #dfx = pd.DataFrame(data=scaledmatrix, columns=columnlabels, index=pcaindex)
    dfx = pd.DataFrame(data=matrix, columns=columnlabels, index=pcaindex)

    # model = pca()
    # The number of components are extracted that cover at least 95% of the explained variance.
    model = pca(n_components=varexp, normalize=True)  # , normalize=True
    results = model.fit_transform(dfx)

    print(model.results['topfeat'])
    outpath = prefix + '.pca.TopFeatures.tsv'
    results['topfeat'].to_csv(outpath, sep='\t')  # pandas dataframe

    # Cumulative explained variance
    print(model.results['explained_var'])
    outpath = prefix + '.pca.ExplainedVariance.tsv'
    results['explained_var'].tofile(outpath)  #numpy.array
    # [0.92461872 0.97768521 0.99478782]

    # Explained variance per PC
    print(model.results['variance_ratio'])
    outpath = prefix + '.pca.VarianceRatio.tsv'
    results['variance_ratio'].tofile(outpath)
    # [0.92461872, 0.05306648, 0.01710261]

    print(model.results['loadings'])
    outpath = prefix + '.pca.Loadings.tsv'
    results['loadings'].to_csv(outpath, sep='\t')

    print(model.results['PC'])
    outpath = prefix + '.pca.PC.tsv'
    results['PC'].to_csv(outpath, sep='\t')

    print(model.results['outliers'])
    outpath = prefix + '.pca.Outliers.tsv'
    results['outliers'].to_csv(outpath, sep='\t')

    # Make plot
    fig, ax = model.plot()
    outpath = prefix + '.pca.VarianceExplained.pdf'
    print(outpath)
    fig.savefig(outpath)  ##### use FIG.savfig, not PLT.savefig... #####
    plt.show()
    plt.close()

    # 2D plot
    fig, ax = model.scatter()
    outpath = prefix + '.pca.scatter.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()
    # 3d Plot
    fig, ax = model.scatter3d()
    outpath = prefix + '.pca.scatter.3D.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    # 2D plot
    fig, ax = model.biplot(n_feat=7, PC=[0, 1])  # n_feat controls number of vectors
    outpath = prefix + '.pca.biplot.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()
    # 3d Plot
    fig, ax = model.biplot3d(n_feat=7, PC=[0, 1, 2])
    outpath = prefix + '.pca.biplot.3D.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    return


def normalize_clusters_for_histogram_sense_antisense(dclusters, clusternames, dgff3, gff3names, nbins=20):
    #  -> Tuple[dict, list]
    import matplotlib.pyplot as plt
    import numpy as np

    # dnorm is a dictionary
    # dnorm key is a float (midpoint of range 0.00-0.05, ...)
    # dnorm value is list [sumS, sumA], sumS is sum of all sense reads for this bin, and same for sumA for antisense
    bins, halfbins = make_histogram_bins(nbins)
    dnorm = {midpoint: [0, 0] for midpoint in halfbins}  # initiate one key per bin for histogram with frequency value 0
    for name in clusternames:
        scaff = name.split('|')[0]
        start = int(name.split('|')[1])
        end = int(name.split('|')[2])
        geneids, geneorientations, genebinary = gff3_gene_overlap(dgff3, gff3names, scaff, start, end)
        # only gather information for clusters overlapping genes # gid is [] if no gene overlap
        for orientation in geneorientations:
            coordinates = dclusters[name][0]  # list of coordinates
            flst = dclusters[name][1]  # list of number of forward reads per coordinate
            rlst = dclusters[name][2]  # list of number of reverse reads per coordinate
            clusterlength = len(coordinates)
            for i in range(clusterlength):  # coordinates are actual scaffold positions, we need 0 ... len(cluster)
                # long clusters will have more entries / bin than small clusters so divide by clusterlength
                f = flst[i] / clusterlength
                r = rlst[i] / clusterlength
                # make cluster positions relative 0-1
                normpos = i / clusterlength
                for count, b in enumerate(bins):  # for each histogram bin, add the f and r values per normposition
                    left, right = b[0], b[1]
                    k = halfbins[count]  # k is name and midpoint of bin (0.050 if bin is 0.0-0.1)
                    if (normpos >= left) and (normpos < right):
                        # f is positive, r is negative. I want sense positive and antisense negative
                        if orientation == '+':  # sense reads in 0 position, antisense in 1 position
                            dnorm[k] = [dnorm[k][0] + f, dnorm[k][1] + r]
                        elif orientation == '-':
                            dnorm[k] = [dnorm[k][0] - r, dnorm[k][1] - f]  # -r bc -(-number) == + and -f bc -(+num) = -
                    # only need to add f and r values to one histogram bin

    return dnorm, halfbins


def make_histogram_bins(nbins, numdecimals=4):
    bins = []
    halfbins = []
    step = 1 / nbins
    halfstep = step / 2
    count = 0
    for i in range(nbins):
        bins.append([round(count, numdecimals), round(count + step, numdecimals)])
        halfbins.append(round(count + halfstep, numdecimals))
        count += step

    return bins, halfbins


def normalize_clusters_for_histogram(dclusters, clusternames, nbins=20):
    # -> Tuple[dict, list]
    import matplotlib.pyplot as plt
    import numpy as np

    # dnorm is a dictionary
    # dnorm key is a float representing the midpoint of a range, ex: 0.025 for 0.00-0.05, ...)
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse
    bins, halfbins = make_histogram_bins(nbins)
    # initiate one key per bin for histogram with frequency value 0 and 0 for f and r
    dnorm = {midpoint: [0, 0] for midpoint in halfbins}
    for name in clusternames:
        coordinates = dclusters[name][0]  # list of coordinates
        flst = dclusters[name][1]  # list of number of forward reads per coordinate
        rlst = dclusters[name][2]  # list of number of reverse reads per coordinate
        clusterlength = len(coordinates)
        for i in range(clusterlength):  # coordinates are actual scaffold positions, we need 0 ... len(cluster)
            # long clusters will have more entries / bin than small clusters so divide by clusterlength
            f = flst[i] / clusterlength
            r = rlst[i] / clusterlength
            # make cluster positions relative 0-1
            normpos = i / clusterlength
            for count, b in enumerate(bins):  # for each histogram bin, add the f and r values per normposition
                left, right = b[0], b[1]
                k = halfbins[count]  # round(start + halfstep, 4) # k is name and midpoint of bin (0.05 if bin 0.0-0.1)
                if (normpos >= left) and (normpos < right):
                    dnorm[k] = [dnorm[k][0] + f, dnorm[k][1] + r]
                # only need to add f and r values to one histogram bin

    return dnorm, halfbins


def plot_per_base_strandedness_all_clusters_concatenated(dclusters, clusternames, prefix, color1='blue', color2='orange', label1='F', label2='R', line_width=1.5):
    import matplotlib.pyplot as plt
    import os

    path, plottitle = os.path.split(prefix)

    position = 0
    x = []
    y1 = []
    y2 = []
    clusterends = []
    for name in clusternames:
        coordinates, sense, antisense = dclusters[name]
        start = position
        end = position + len(coordinates)
        for i in range(start, end):
            x.append(i)
        for s, a in zip(sense, antisense):
            y1.append(s)
            y2.append(a*-1)
        position = end
        clusterends.append(end)

    plt.figure(figsize=(15, 7.5))
    plt.plot(x, y1, color=color1, label=label1, linewidth=line_width)  # 'F' or 'S'  # 0.5
    plt.plot(x, y2, color=color2, label=label2, linewidth=line_width)  # # 'R' or 'A' # 0.5
    plt.grid()
    plt.legend(loc='upper left')
    plt.title(plottitle)

    outpath = prefix + '.perbase.strandedness.AllClusters.Concatenated.plots.pdf'
    print(outpath)
    plt.savefig(outpath)
    plt.show()
    plt.close()

    return


def plot_per_base_strandedness(dclusters, clusternames, prefix, dbed, bednames, numplots=4, color1='blue', color2='orange', label1='F', label2='R'):
    import matplotlib.pyplot as plt

    for i in range(0, len(clusternames), numplots):
        numsubplots = int(numplots / 2)
        fig, ax = plt.subplots(numsubplots, numsubplots, figsize=(10, 10))  # sharex=False
        a, b = 0, 0
        for j in range(i, i + numplots):
            if j == len(clusternames):
                break  # break to avoid out of range clusternames[j]
            if j == i:
                pass
            elif j == i + numplots:
                b += 1
            elif j % 2 == 1:
                b += 1
            elif j % 2 == 0:
                a += 1
                b -= 1
            name = clusternames[j]
            scaff = name.split('|')[0]
            start = name.split('|')[1]
            end = name.split('|')[2]
            # srcids is list of all src cluster ids that overlap window
            srcids = src_bed_overlap(dbed, bednames, scaff, int(start), int(end))
            ids = ','.join(srcids)
            if ids == '':
                ids = 'NewCluster'
            x = dclusters[name][0]  # genomic positions, x axis coordinates
            y = dclusters[name][1]  # number of forward reads per x axis coordinate
            z = dclusters[name][2]  # number of reverse reads per x axis coordinate
            z = [-1*number for number in z]  # make negative
            ax[a, b].plot(x, y, color=color1, label=label1)
            ax[a, b].plot(x, z, color=color2, label=label2)
            ax[a, b].set_title(ids + ': ' + name[len('scaffold51_'):])
            ax[a, b].set_xlim([int(start), int(end)])
            ax[a, b].legend(loc='upper left')
            for k, label in enumerate(ax[a, b].get_xticklabels()):
                if k % 2 == 1:
                    label.set_visible(True)
                elif k % 2 == 0:
                    label.set_visible(False)
        outpath = prefix + '.perbase.strandedness.plots' + str(i) + '.pdf'
        print(outpath)
        plt.savefig(outpath)
        plt.show()
        plt.close()

    return


def normalize_per_base_depth_by_allreads(dclusters, allreadcount) -> dict:
    # dclusters key is 'scaff|start|end'
    # value is [[coordinates], [forwardreadcounts], [reversereadcounts]], lists contain int()
    dclustersnorm = {}
    for k, v in dclusters.items():
        p, f, r = v
        newf = [(i / allreadcount) * 1_000_000 for i in f]
        newr = [(i / allreadcount) * 1_000_000 for i in r]
        dclustersnorm[k] = [p, newf, newr]

    return dclustersnorm


def per_base_strandedness_orient_by_gene_limit(bamfile, strandedfilename, ids):
    # ids list of ids from limiting bed file  # class I SRCs
    # windows = [scaffold, start, end]
    # -> Tuple[dict, list]
    import pysam

    with open(strandedfilename) as FILE:
        header = FILE.readline()
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        srcid = window[0]
        scaff = window[1]
        start = int(window[2])
        end = int(window[3])
        polarity = window[-13]  # SRC_VS_gene_polarity 'A', 'S', or '-'. not desired 'X' or -,-, or A,S
        gene_orientation = window[-14]  # gene_orientations  # + or -. not desired NA, +,_ , +,+ , -,-
        name = '|'.join([scaff, str(start), str(end), polarity])
        if srcid in ids:

            f, r, coordinates = per_base_strandedness_counts(samfile, scaff, start, end)

            if gene_orientation == '+':  # if mRNA is stranded in direction of genome '+'
                dclusters[name] = [coordinates, f, r]
                clusternames.append(name)
            elif gene_orientation == '-':  # if mRNA is stranded in opposite direction of genome '-'
                # flip & swap per base depth values to match mRNA orientation
                dclusters[name] = [coordinates, r[::-1], f[::-1]]
                clusternames.append(name)

    # dclusters has all sense or antisense clusters oriented in direction of mRNA
    # all mRNA are oriented to be '+' in direction of genome.
    # here f is sense, r is antisense
    # only append cluster name if sense or antisense to mRNA

    return dclusters, clusternames


def per_base_strandedness_orient_by_mrna_limit(bamfile, strandedfilename, ids):
    # ids list of ids from limiting bed file
    # windows = [scaffold, start, end]
    import pysam

    with open(strandedfilename) as FILE:
        header = FILE.readline()
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        srcid = window[0]
        scaff = window[1]
        start = int(window[2])
        end = int(window[3])
        mrnafor = float(window[-5])  # mRNA forward percent
        polarity = window[-4]  # SRC_vs_mRNA_polarity 'A', 'S', or 'N'
        name = '|'.join([scaff, str(start), str(end), polarity])
        if srcid in ids:

            f, r, coordinates = per_base_strandedness_counts(samfile, scaff, start, end)

            if mrnafor >= 75.0:  # if mRNA is stranded in direction of genome '+'
                dclusters[name] = [coordinates, f, r]
                clusternames.append(name)

            elif mrnafor <= 25.0:  # if mRNA is stranded in opposite direction of genome '-'
                # flip & swap per base depth values to match mRNA orientation
                dclusters[name] = [coordinates, r[::-1], f[::-1]]
                clusternames.append(name)
            #dclusters[name] = [coordinates, f, r]

    # dclusters has all sense or antisense clusters oriented in direction of mRNA
    # all mRNA are oriented to be '+' in direction of genome.
    # here f is sense, r is antisense
    # only append cluster name if sense or antisense to mRNA

    return dclusters, clusternames


def per_base_strandedness_orient_by_gene(bamfile, strandedfilename):
    # -> Tuple[dict, list]
    import pysam

    with open(strandedfilename) as FILE:
        header = FILE.readline()
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        scaff = window[1]
        start = int(window[2])
        end = int(window[3])
        gc_reads = float(window[5])
        polarity = window[-13]  # SRC_VS_gene_polarity 'A', 'S', or '-'. not desired 'X' or -,-, or A,S
        gene_orientation = window[-14]  # gene_orientations  # + or -. not desired NA, +,_ , +,+ , -,-
        name = '|'.join([scaff, str(start), str(end), polarity])
        if (polarity == 'S') or (polarity == 'A') or (polarity == '-'):
            clusternames.append(name)

            f, r, coordinates = per_base_strandedness_counts(samfile, scaff, start, end)

            if gene_orientation == '+':  # if mRNA is stranded in direction of genome '+'
                dclusters[name] = [coordinates, f, r]
            elif gene_orientation == '-':  # if mRNA is stranded in opposite direction of genome '-'
                # flip & swap per base depth values to match mRNA orientation
                dclusters[name] = [coordinates, r[::-1], f[::-1]]

    # dclusters has all sense or antisense clusters oriented in direction of mRNA
    # all mRNA are oriented to be '+' in direction of genome.
    # here f is sense, r is antisense
    # only append cluster name if sense or antisense to mRNA

    return dclusters, clusternames


def per_base_strandedness_orient_by_mrna(bamfile, strandedfilename):
    import pysam

    with open(strandedfilename) as FILE:
        header = FILE.readline()
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        scaff = window[1]
        start = int(window[2])
        end = int(window[3])
        gc_reads = float(window[5])
        mrnafor = float(window[-5])  # mRNA forward percent
        polarity = window[-4]  # SRC_vs_mRNA_polarity 'A', 'S', or 'N'
        name = '|'.join([scaff, str(start), str(end), polarity])
        if (polarity == 'S') or (polarity == 'A'):
            clusternames.append(name)

            f, r, coordinates = per_base_strandedness_counts(samfile, scaff, start, end)

            if mrnafor >= 75.0:  # if mRNA is stranded in direction of genome '+'
                dclusters[name] = [coordinates, f, r]
            elif mrnafor <= 25.0:  # if mRNA is stranded in opposite direction of genome '-'
                # flip & swap per base depth values to match mRNA orientation
                dclusters[name] = [coordinates, r[::-1], f[::-1]]
            #dclusters[name] = [coordinates, f, r]

    # dclusters has all sense or antisense clusters oriented in direction of mRNA
    # all mRNA are oriented to be '+' in direction of genome.
    # here f is sense, r is antisense
    # only append cluster name if sense or antisense to mRNA

    return dclusters, clusternames


def per_base_strandedness_counts(samfile, scaff, start, end):
    import pysam

    # samfile = pysam.AlignmentFile(bamfile, "rb")  # it needs to already be opened with pysam
    coordinates, f, r = [], [], []  # f = sum of sense reads per base, r = sum of antisense reads per base
    # positions = [i for i in range(start, end)]
    # for pos in range(start, end):
    for pileupcolumn in samfile.pileup(scaff, start, end, stepper='nofilter', max_depth=100_000):
        # pileup gives depth values outside region
        # if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
        if (pileupcolumn.pos <= end) and (pileupcolumn.pos >= start):
            forward, reverse = 0, 0
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.alignment.is_reverse == True:
                        reverse += 1
                    elif pileupread.alignment.is_reverse == False:
                        forward += 1
            coordinates.append(pileupcolumn.pos)
            f.append(forward)
            r.append(reverse)

    return f, r, coordinates


def per_base_strandedness(bamfile, windows):
    # windows = [SRCID, scaffold, start, end]
    # windowfile
    import pysam

    # with open(windowfile) as FILE:
    #     windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        name = '|'.join([scaff, window[1], window[2]])
        clusternames.append(name)
        f, r, coordinates = per_base_strandedness_counts(samfile, scaff, start, end)
        dclusters[name] = [coordinates, f, r]

    return dclusters, clusternames


def update_src_bed_file_with_new_srcs(srcbedfile, newclusters):
    with open(srcbedfile, 'r') as FILE:
        srcs = [line.strip() for line in FILE if line.strip() != ''] + newclusters
    with open(srcbedfile, 'w') as OUT:
        OUT.write('\n'.join(srcs) + '\n')
    return


def count_new_srcs_compared_to_paper(srcfilefrompaper, strandwindows):
    with open(srcfilefrompaper, 'r') as FILE:
        header = FILE.readline()
        srcids = [line.strip().split('\t')[3] for line in FILE if line.strip() != '']

    identifiedsrcids = [row.split('\t')[0] for row in strandwindows[1:]]  # skip first row header # first entry in each row is SRCID
    # remove scrids from identifiednsrcids
    numbernewscridsrelativetopaper = len(set(identifiedsrcids) - set(srcids))
    print(', '.join(list(set(identifiedsrcids) - set(srcids))))

    print(f'Number of new SRC IDs relative to Karunanithi et al. 2019: {numbernewscridsrelativetopaper}')
    return


def window_read_count_stranded_mRNA(samfile, scaff, start, end, coveredbymRNA):
    # samfile was already opened with (samfile = pysam.AlignmentFile(bamfile, "rb")
    # flag is integer value
    # for truseq mrna stranded kit!
    countplus = 0  #f2r1
    countminus = 0  #f1r2
    plusbinary = 0
    minusbinary = 0
    for read in samfile.fetch(scaff, start, end):
        line = read.tostring()
        flag = int(line.strip().split('\t')[1])
        # only count second reads in pairs # 1 count per fragment
        # read paired, read in proper pair, not unmapped, mate not unmapped, read not on reverse strand, mate on reverse strand, not first in pair, second in pair
        if f'{flag:012b}'[-8:] == '10100011':  # f2r1 # coding strand is + (mRNA same orientation as genome)
            countplus += 1
        # read paired, read in proper pair, not unmapped, mate not unmapped, read on reverse strand, mate not on reverse strand, not first in pair, second in pair
        elif f'{flag:012b}'[-8:] == '10010011':  # f1r2 # coding strand is - (mRNA opposite orientation as genome)
            countminus += 1
    total = countplus + countminus
    normplus = (countplus / (end - start + 1)) * 1000  # 1000 is a scaling factor to make values more readable
    normminus = (countminus / (end - start + 1)) * 1000  # 1000 is a scaling factor to make values more readable
    normbysizeplus = round(normplus, 2)  # normalizes mRNA count by cluster size
    normbysizeminus = round(normminus, 2)  # normalizes mRNA count by cluster size
    normtotal = normbysizeplus + normbysizeminus

    # strandplus = round((normbysizeplus / normtotal) * 100, 2)
    if (total > 0) and (coveredbymRNA == 1):  # (total > 0) and (normtotal >= 5)
        strandplus = round((countplus / total) * 100, 2)  # percent of mRNA reads in + orientation
        if strandplus >= 75.0:
            plusbinary = 1
        elif strandplus <= 25.0:
            minusbinary = 1
    elif (total > 0):  # (total > 0) and (normtotal >= 5)
        strandplus = round((countplus / total) * 100, 2)  # percent of mRNA reads in + orientation
        if strandplus >= 50.0:
            strandplus = 51.0
        elif strandplus < 50.0:
            strandplus = 49.0
    else:
        strandplus = 50.0

    return countplus, countminus, normbysizeplus, normbysizeminus, strandplus, plusbinary, minusbinary


def window_read_count_rpkm(samfile, scaff, start, end, librarysize):
    import pysam
    # samfile was already opened with (samfile = pysam.AlignmentFile(bamfile, "rb")
    count = 0
    for read in samfile.fetch(scaff, start, end):
        count += 1
    norm = (count / (end - start + 1)) * 1_000  # 1000 is a scaling factor to make values more readable
    countrpkm = round((norm / librarysize) * 1_000_000, 2)  # 1000000 is scaling factor to make rpkm
    normbylength = round(norm, 2)  # normalizes mRNA count by cluster size
    if normbylength >= 100:
        coveredbymRNA = 1
    else:
        coveredbymRNA = 0

    return count, normbylength, countrpkm, coveredbymRNA


def window_read_count_forward_reverse(samfile, scaff, start, end):
    import pysam
    forward, reverse = 0, 0
    for read in samfile.fetch(scaff, start, end):
        if read.is_reverse == True:
            reverse += 1
        elif read.is_reverse == False:
            forward += 1
    totalreads = forward + reverse

    return forward, reverse, totalreads


def define_sRNA_gene_polarity(strandfor, strandrev, geneids, geneorientations):
    polarity = []
    if (strandfor >= 75.0) or (strandrev >= 75.0):
        strandbinary = 1  # 0 if cluster is not stranded (<75% reads of one orientation) 1 if stranded
    else:
        strandbinary = 0
    # 0 if cluster reads are sense to all coevered genes, or cluster reads not stranded
    # 1 if cluster reads are antisense to at least one covered gene
    # if multiple genes cover cluster, value is biased to 1, if one gene is antisense return 1
    antisensebinary = 0
    gids = ','.join(geneids)
    orientations = ','.join(geneorientations)
    if gids == '':
        gids = 'NA'  # no gene overlap
        orientations = 'NA'
        polarity.append('X')  # polarity will be S, A, N, -, or X
    else:
        for gorientation in geneorientations:
            if strandfor >= 75.0:
                if gorientation == '+':
                    polarity.append('S')  # S == sense
                elif gorientation == '-':
                    polarity.append('A')  # A == antisense
                    antisensebinary = 1
                elif gorientation == '.':
                    polarity.append('N')  # N == no gene orientation information
            elif strandrev >= 75.0:
                if gorientation == '+':
                    polarity.append('A')  # S == sense
                    antisensebinary = 1
                elif gorientation == '-':
                    polarity.append('S')  # A == antisense
                elif gorientation == '.':
                    polarity.append('N')  # N == no gene orientation information
            else:  # if if strandfor and strandrev are both < 75.0 (no majority of reads pointing in same direction)
                polarity.append('-')  # - == no strong strand bias for cluster (can't say sense or antisense)
    pols = ','.join(polarity)

    return strandbinary, antisensebinary, orientations, gids, pols


def gff3_gene_overlap(dgff3, gff3names, scaff, start, end, coveragecutoff=75.0):
    geneids = []
    orientations = []
    genebinary = 0  # 0 if overlaps no genes, 1 if overlaps at least one gene
    for k in gff3names:
        genescaff = dgff3[k][0]
        genestart = int(dgff3[k][3])
        geneend = int(dgff3[k][4])
        geneorientation = dgff3[k][6]
        id = dgff3[k][8].split(';')[0][3:]  # isolate gene id, from metadata column
        if (scaff == genescaff) and (start <= geneend) and (end >= genestart):
            # below calculates amount of overlap, if they don't overlap, returns 0
            # len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1))
            overlap = len(range(max(start, genestart), min(end, geneend) + 1))
            percentoverlap = (overlap / (end - start + 1)) * 100
            if percentoverlap >= coveragecutoff:  # if gene covers more than 25% of cluster record it
                geneids.append(id)
                orientations.append(geneorientation)
                genebinary = 1

    return geneids, orientations, genebinary


def src_bed_overlap(dbed, bednames, scaff, start, end):
    srcids = []
    for k in bednames:
        bedscaff = dbed[k][0]
        bedstart = int(dbed[k][1])
        bedend = int(dbed[k][2])
        id = dbed[k][3]
        # if (scaff == bedscaff) and (start <= bedend) and (end >= bedstart):
        if (scaff == bedscaff):
            # below calculates amount of overlap, if they don't overlap, returns 0
            # len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1))
            percent_overlap = (len(range(max(start, bedstart), min(end, bedend) + 1)) / (end - start + 1)) * 100
            if percent_overlap >= 25:
                srcids.append(id)

    return srcids


def fasta_seq_gc_percent(d, scaff, start, end):
    GCcount = 0
    seq = d[scaff][start:end+1]
    for base in seq:
        if (base == 'G') or (base == 'C'):
            GCcount += 1
    refGC = round((GCcount / len(seq)) * 100, 2)

    return refGC, seq


def count_all_reads_in_bam(bamfile):
    import pysam
    # idxstats
    # output is TAB-delimited with each line consisting of:
    # reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments
    # print reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(filename) ])
    stats = pysam.idxstats(bamfile)
    statslines = stats.strip().split('\n')  # strip removes trailing \n
    allreads = sum([int(l.split('\t')[2]) + int(l.split('\t')[3]) for l in statslines])
    mappedreads = sum([int(l.split('\t')[2]) for l in statslines])
    unmappedreads = sum([int(l.split('\t')[3]) for l in statslines])

    return allreads, mappedreads, unmappedreads


def calculate_read_gc_percent(samfile, scaff, winstart, winend):
    # samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    import pysam
    import numpy as np

    coverage = samfile.count_coverage(scaff, winstart, winend, read_callback='nofilter')
    numGCs = np.sum(coverage[1] + coverage[2])  # total number of GC bases for all aligned reads in window
    numbases = np.sum(coverage)
    if not numbases:  # if zero
        numbases = 1
    GCpercent = round((numGCs / numbases) * 100, 2)  # GC percentage of all aligned reads

    return GCpercent


def calculate_average_depth_window(samfile, scaff, winstart, winend):
    # samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    import pysam
    import numpy as np

    coverage = samfile.count_coverage(scaff, winstart, winend, read_callback='nofilter')
    depth = list(np.sum(coverage, axis=0))  # sum four arrays vertically with np.sum axis=0
    winlength = winend - winstart + 1
    average_depth = mean_list(depth, denominator=winlength)

    return average_depth


def strandedness(bamfile, mRNAfile, windows, dbed, bednames, dgff3, gff3names, dgenome, max_id, srnaallreads, mrnaallreads):
    # function suffers from scope-creep
    # windowfile
    import pysam

    # with open(windowfile) as FILE:
    #     windows = [line.strip().split('\t') for line in FILE]

    newwindows = []
    pcatable = []
    pcarownames = []
    pcaindex = []
    newclusters = []
    countnewcluster = 0
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    rnafile = pysam.AlignmentFile(mRNAfile, "rb")
    # srnaallreads, srnamappedreads, srnaunmappedreads = count_all_reads_in_bam(bamfile)
    # mrnaallreads, mrnamappedreads, mrnaunmappedreads = count_all_reads_in_bam(mRNAfile)
    for window in windows:
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        # averagedepth = round(float(window[3]), 2)
        # readsgc = round(float(window[4]), 2)
        averagedepth = round(calculate_average_depth_window(samfile, scaff, start, end), 2)
        readsgc = calculate_read_gc_percent(samfile, scaff, start, end)
        window.append(str(averagedepth))
        window.append(str(readsgc))

        length = end - start + 1
        refseqGCpercent, seq = fasta_seq_gc_percent(dgenome, scaff, start, end)

        forward, reverse, totalreads = window_read_count_forward_reverse(samfile, scaff, start, end)

        sRNA_RPKM = round((totalreads / (length * srnaallreads)) * 1_000 * 1_000_000, 2)  # 1000 is a scaling factor to make values more readable
        strandfor = round((forward / totalreads) * 1_00, 2)  # percent of reads in forward direction
        strandrev = round((reverse / totalreads) * 1_00, 2)  # percent of reads in reverse direction
        # srcids is list of all src cluster ids that overlap window
        srcids = src_bed_overlap(dbed, bednames, scaff, start, end)
        ids = ','.join(srcids)
        if ids == '':  # if not ids
            countnewcluster += 1
            ids = f'C{countnewcluster + max_id}'  # id format c1, c2, c2, ... , c421 ...
            newclusters.append('\t'.join([scaff, str(start), str(end), f'C{countnewcluster + max_id}']))

        # determine if sRNA overlap gene
        geneids, geneorientations, genebinary = gff3_gene_overlap(dgff3, gff3names, scaff, start, end)
        # determine sRNA polarity to gene
        strandbinary, antisensebinary, orientations, gids, pols = define_sRNA_gene_polarity(strandfor, strandrev, geneids, geneorientations)

        # get number of mRNA reads overlapping window
        mRNAcount, countnormbylength, mRNA_RPKM, coveredbymRNA = window_read_count_rpkm(rnafile, scaff, start, end, mrnaallreads)  # raw count
        countplus, countminus, normbysizeplus, normbysizeminus, strandplus, plusbinary, minusbinary = window_read_count_stranded_mRNA(rnafile, scaff, start, end, coveredbymRNA)  # raw count

        countnormplusminus = normbysizeplus + normbysizeminus
        mrnapolarity = ''  # see comments below
        # strandplus is mRNA percent in forward direction
        # strandfor is sRNA percent in forward direction
        # N if mRNA or sRNA over SRC are not stranded
        # S if sRNA from SRC are sense to mRNA reads over SRC
        # A if sRNA from SRC are antisense to mRNA reads over SRC
        if (strandplus >= 75.0) and (strandfor >= 75.0):  # used 5 and not 10 because i only count reverse reads, therefore counts will be ~half of all counts.
            mrnapolarity = 'S'
        elif (strandplus <= 25.0) and (strandfor <= 25.0):
            mrnapolarity = 'S'
        elif (strandplus >= 75.0) and (strandfor <= 25.0):
            mrnapolarity = 'A'
        elif (strandplus <= 25.0) and (strandfor >= 75.0):
            mrnapolarity = 'A'
        else:
            mrnapolarity = 'N'

        newwindows.append('\t'.join([ids] + window + [str(refseqGCpercent), str(length), str(sRNA_RPKM), str(forward), str(reverse), str(strandfor), str(strandrev), str(strandbinary), gids, orientations, pols, str(antisensebinary), str(mRNAcount), str(countnormbylength), str(mRNA_RPKM), str(countplus), str(countminus), str(countnormplusminus), str(strandplus), mrnapolarity, str(plusbinary), str(minusbinary), seq]))
        pcatable.append([str(sRNA_RPKM), str(readsgc), str(refseqGCpercent), str(length), str(strandfor), str(mRNA_RPKM), str(strandplus)])  # str(genebinary), str(antisensebinary),
        pcarownames.append(ids)
        pcaindex.append(str(coveredbymRNA))
    header = '\t'.join(['SRC_ID', 'scaffold', 'start', 'end', 'average_depth', 'GCPercent_reads', 'GCPercent_reference', 'length', 'SRC_RPKM', 'forward_reads', 'reverse_reads', 'forward_percent', 'reverse_percent', 'stranded',  'geneIDs', 'gene_orientations', 'SRC_VS_gene_polarity', 'antisense', 'mRNA_count', 'mRNANormBylength', 'mRNA_RPKM', 'mRNA_plus_count', 'mRNA_minus_count', 'mRNA_stranded_count_normbySRClength', 'mRNA_forward_percent', 'SRC_VS_mRNA_polarity', 'mRNA_stranded_forward', 'mRNA_stranded_reverse', 'sequence'])
    pcacolumnnames = ['SRC_RPKM', 'ReadGC', 'ReferenceGC', 'length', 'SRC_forward_percent', 'mRNA_RPKM', 'mRNA_forward_percent']  # 'covered_by_gene', 'antisense',
    samfile.close()

    outwindows = [header] + newwindows
    # pcaoutput = [pcaheader] + pcatable

    # x = [1, 2, 3]
    # y = [[4, 5, 6], [7, 8, 9]]
    # x + y = [1, 2, 3, [4, 5, 6], [7, 8, 9]]
    # [x] + y = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

    print(f'Number of new clusters identified: {countnewcluster}')
    print(f'Number of total clusters identified: {len(newwindows)}')

    return outwindows, pcatable, pcacolumnnames, pcarownames, pcaindex, newclusters


def refine_window_coordinates(sortwindows, bamfile) -> list:
    import pysam
    import numpy as np

    refinedwindows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    for win in sortwindows:
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        # array.arrays A, C, G, T. bases, per coordinate
        coverage = samfile.count_coverage(scaff, winstart, winend, read_callback='nofilter')
        numbases = np.sum(coverage)
        depth = list(np.sum(coverage, axis=0))  # sum four arrays vertically with np.sum axis=0
        # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
        numGCs = np.sum(coverage[1] + coverage[2])  # total number of GC bases for all aligned reads in window
        countstart, countend = 0, 0
        while depth[0] == 0:
            countstart += 1
            del depth[0]
        while depth[-1] == 0:
            countend += 1
            del depth[-1]
        refinedstart = winstart + countstart
        refinedend = winend - countend
        refinedwinlength = refinedend - refinedstart + 1
        average_depth = mean_list(depth, denominator=refinedwinlength)
        GCpercent = round((numGCs / numbases) * 100, 2)  # GC percentage of all aligned reads
        refinedwindows.append(f'{scaff}\t{refinedstart}\t{refinedend}\t{average_depth}\t{GCpercent}')

    samfile.close()

    #refinedwindows = ['\t'.join(lst) for lst in refinedwindows]
    return refinedwindows


def write_table(table, outpath):
    # tale is list of lists
    print('Writing table to: %s' % outpath)
    output = ['\t'.join(l) for l in table]
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(output))


def write_to_bed(windows, outpath):
    print('Writing windows to: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(windows))


def collapse_windows(sortwindows):
    # this function is janky
    # sortwindows is a naturally sorted list of strings
    # the strings inside the list are formatted as scaffold\tstartposition\tendposition
    collapsedsortwindows = []
    collapsewindepth = []
    collapseGC = []
    overlapgate = 1  # record start position of first window in overlap
    for count, win in enumerate(sortwindows):
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        windepth = float(win.split('\t')[3])
        GCpercent = float(win.split('\t')[4])
        if count == 0:
            pass
        elif (scaff == prevscaff) and (winstart <= prevwinend) and (prevwinstart <= winend):
            # some type of overlap # also equivalent to (i think) a[1] > b[0] and a[0] < b[1]
            # determine if there is overlap between previous window and current window, then combine them
            if overlapgate == 1:
                collapsewindow = '%s\t%s' % (scaff, prevwinstart)  # will add end coordinate when i find end of overlap
                overlapgate = 0
                collapsewindepth.append(prevwindepth)
                collapsewindepth.append(windepth)
                collapseGC.append(prevwinGC)
                collapseGC.append(GCpercent)
        elif overlapgate == 0:
            # if it makes it to this point, there is no overlap with previous window
            # but if overlapgate == 0 then there was overlap between at least two windows and we combine coordinates
            collapsedsortwindows.append('%s\t%s\t%.2f\t%.2f' % (collapsewindow, prevwinend, sum(collapsewindepth) / len(collapsewindepth), sum(collapseGC) / len(collapseGC)))
            collapsewindepth = []  # reset for next series of overlapping windows
            collapseGC = []
            overlapgate = 1  # reset overlapgate to allow for first start coordinates to be recorded
        elif overlapgate == 1:
            # if no overlap with previous, and overlapgate == 1 then no windows were overlapping
            # check if this window WILL overlap with the NEXT window and add it to collapsedsortwindows if NO
            if win != sortwindows[-1]:  # need to have this to avoid list indexing error
                if (scaff == sortwindows[count+1].split('\t')[0]) and (winstart <= int(sortwindows[count+1].split('\t')[2])) and (int(sortwindows[count+1].split('\t')[1]) <= winend):
                    pass
                else:
                    # no overlap with NEXT window
                    collapsedsortwindows.append(win)
                    collapsewindepth = []
                    collapseGC = []
            else:  # if it is the last window here output the window
                collapsedsortwindows.append(win)
                collapsewindepth = []
                collapseGC = []

        prevscaff = win.split('\t')[0]
        prevwinstart = int(win.split('\t')[1])
        prevwinend = int(win.split('\t')[2])
        prevwindepth = float(win.split('\t')[3])
        prevwinGC = float(win.split('\t')[4])

    return collapsedsortwindows


def mean_list(lst, denominator):
    mean = sum(lst) / denominator
    return mean


def mean_with_denominator(value, denominator):
    # mean = sum(lst) / denominator
    mean = value / denominator
    return mean


def calculate_depth_with_pysam_sliding_windows(bamfile, dscaff_max, step=25, winsize=100, cutoff=100):
    import pysam
    import numpy as np

    windows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  #pysam.AlignmentFile("ex1.bam", "rb")
    for scaff, maxcoord in dscaff_max.items():
        print(scaff, end=', ')
        for i in range(0, maxcoord, step):
            ## tempdepth = []  # temppos = []
            ## .pileup skips bases when depth = 0  # need to figure out mpileup -a option in pysam
            ## tempdepth = [pileupcolumn.nsegments for pileupcolumn in samfile.pileup(scaff, i, i + winsize)]
            #tempdepth = []
            #for pos in range(i, i + winsize):
            #    depth = 0
            #    for pileupcolumn in samfile.pileup(scaff, pos, pos+1):  # pileup gives depth values outside region
            #        if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
            #            #for pileupread in pileupcolumn.pileups:
            #            #    if not pileupread.is_del and not pileupread.is_refskip:
            #            #        depth += 1
            #            #tempdepth.append(depth)
            #            ##tempdepth.append(pileupcolumn.nsegments)  # pileupcolumn.nsegments # n depricated
            coverage = samfile.count_coverage(scaff, i, i + winsize, read_callback='nofilter')  # array.arrays A, C, G, T. bases, per coordinate
            numbases = np.sum(coverage)
            # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
            numGCs = np.sum(coverage[1] + coverage[2])  # number of GC bases in window

            if numbases >= 10:
                average_depth = mean_with_denominator(numbases, denominator=winsize)
                GCpercent = round((numGCs / numbases) * 100, 2)
            else:
                average_depth = 0
                GCpercent = 0
            if average_depth >= cutoff:
                # print(f'{average_depth}, {scaff}, {i+1}, {i+winsize}, {tempdepth}')
                # print(str(average_depth), end=', ')
                # scaff, start, end, averagedepth
                windows.append([scaff, str(i + 1), str(i + winsize), str(average_depth), str(GCpercent)])
                for j in range(i + 1, i + winsize):
                    if j <= maxcoord:
                        positions.add('_'.join([scaff, str(j)]))
    samfile.close()

    return windows, positions


def read_gff3(f):
    # assumes gff3 is trimmed to 'gene' lines only (no exon/CDS/whatever else), one line per gene
    dgff3 = {}
    gff3names = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] != '#':
                s = line.strip().split('\t')
                dgff3['|'.join([s[0], s[3], s[4]])] = line.strip().split('\t')
                gff3names.append('|'.join([s[0], s[3], s[4]]))

    return dgff3, gff3names


def read_bed(f):
    dbed = {}
    bednames = []
    bedids = []
    with open(f, 'r') as FILE:
        header = FILE.readline()
        for line in FILE:
            dbed['|'.join(line.strip().split('\t')[:3])] = line.strip().split('\t')
            bednames.append('|'.join(line.strip().split('\t')[:3]))
            bedids.append(line.strip().split('\t')[3])

    return dbed, bednames, bedids


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:-len('_with_IES')]
                    namelist.append(name)
                    d[name] = []
                else:
                    name = line.strip()[1:]
                    namelist.append(name)
                    d[name] = []
            elif line.strip() != '':  # else: # trying to prevent some crap happening on the last line
                d[name].append(line.strip())
    for name in namelist:
        d[name] = ''.join(d[name])  # join list of partial sequences. Useful if interleaved fasta

    return d, namelist


def length_of_fasta_sequences(genomefile):
    import os

    print('Counting lengths of all scaffolds')
    path, f = os.path.split(genomefile)
    dgenome, names = read_fasta_as_dict(genomefile)
    d = {k: len(v) for k, v in dgenome.items()}

    return d, dgenome, names


def max_coords_per_scaffold(d):
    dmax_coords = {}
    dscaff_max = {}
    for n in list(d.keys()):
        # collect all coordinates in list, for each scaffold
        dmax_coords.setdefault('_'.join(n.split('_')[:-1]), []).append(int(n.split('_')[-1]))
    for scaff in list(dmax_coords.keys()):
        # get maximum coordinate per scaffold
        dscaff_max[scaff] = max(dmax_coords[scaff])
    return dscaff_max


def plot_pca_and_pairwise_correlation(pcatable, pcacolumnnames, pcaindex, prefix, varexp):
    # dfx = pd.DataFrame(data=pcatable, columns=pcacolumnnames, index=pcaindex)
    # print(dfx.head())
    make_pca(pcatable, pcacolumnnames, pcaindex, prefix, varexp)
    # pairdf = pd.DataFrame(data=pcatable, columns=pcacolumnnames, index=pcarownames)
    make_sns_pair_plot(pcatable, pcacolumnnames, pcaindex, prefix)

    return


def plot_reads_oriented_by_genome(dclusters, clusternames, srnaallreads, dbed, bednames, numplots, prefix):
    dclustersnorm = normalize_per_base_depth_by_allreads(dclusters, srnaallreads)
    plot_per_base_strandedness(dclustersnorm, clusternames, prefix, dbed, bednames, numplots, color1='darkcyan',
                               color2='red', label1='F', label2='R')
    ####### stranded plot
    temp = prefix + '.ForwardReverseFingerprint'
    plot_per_base_strandedness_all_clusters_concatenated(dclustersnorm, clusternames, temp, color1='darkcyan',
                                                         color2='red', label1='F', label2='R', line_width=0.5)
    dnorm, halfbins = normalize_clusters_for_histogram(dclustersnorm, clusternames, nbins=50)
    temp = prefix + '.Forward.Reverse'
    plot_one_normalized_per_base_strandedness(dnorm, halfbins, temp)
    temp = prefix + '.ForwardMinusReverse'
    plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, temp, color1='darkcyan', label1='F - R')

    return


def plot_reads_oriented_by_mRNA(dclusterspolar, clusternamespolar, srnaallreads, dbed, bednames, numplots, prefix):
    dclusterspolarnorm = normalize_per_base_depth_by_allreads(dclusterspolar, srnaallreads)

    temp = prefix + '.OrientedByForwardmRNA'
    plot_per_base_strandedness(dclusterspolarnorm, clusternamespolar, temp, dbed, bednames, numplots, color1='orange',
                               color2='blue', label1='S', label2='A')
    temp = prefix + '.OrientedByForwardmRNA'
    plot_per_base_strandedness_all_clusters_concatenated(dclusterspolarnorm, clusternamespolar, temp, color1='orange',
                                                         color2='blue', label1='S', label2='A', line_width=0.5)
    dnorm, halfbins = normalize_clusters_for_histogram(dclusterspolarnorm, clusternamespolar, nbins=50)
    temp = prefix + '.Sense.Antisense'
    plot_one_normalized_per_base_strandedness(dnorm, halfbins, temp, color1='orange', color2='blue', label1='S',
                                              label2='A')
    temp = prefix + '.SenseMinusAntisense'
    plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, temp, color1='blue', label1='S - A')

    return


def plot_reads_oriented_by_gene(dclustersgene, clusternamesgene, srnaallreads, dbed, bednames, numplots, prefix):
    dclustersgenenorm = normalize_per_base_depth_by_allreads(dclustersgene, srnaallreads)

    temp = prefix + '.OrientedByForwardGene'
    plot_per_base_strandedness(dclustersgenenorm, clusternamesgene, temp, dbed, bednames, numplots, color1='orange',
                               color2='blue', label1='S', label2='A')
    temp = prefix + '.OrientedByForwardGene'
    plot_per_base_strandedness_all_clusters_concatenated(dclustersgenenorm, clusternamesgene, temp, color1='orange',
                                                         color2='blue', label1='S', label2='A', line_width=1.0)
    dnorm, halfbins = normalize_clusters_for_histogram(dclustersgenenorm, clusternamesgene, nbins=50)
    temp = prefix + '.Sense.Antisense.Gene'
    plot_one_normalized_per_base_strandedness(dnorm, halfbins, temp, color1='orange', color2='blue', label1='S',
                                              label2='A')
    temp = prefix + '.SenseMinusAntisense.Gene'
    plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, temp, color1='blue', label1='S - A')

    return


def main_limit_clusters(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step=25, winsize=100, cutoff=100,
                        numplots=4, varexp=0.9, libsize=0, limitbedall='', limitbed1='', limitbed2=''):
    import os
    from natsort import natsorted, ns
    import pandas as pd
    ''' limitbedall is tab separated bed file, four columns, scaff\tstart\tend\tSRCID'''
    ''' limitbed1 is tab separated bed file, four columns, subset of limitbedall, Class1, scaff\tstart\tend\tSRCID'''
    ''' limitbed2 is tab separated bed file, four columns, subset of limitbedall, Class2, scaff\tstart\tend\tSRCID'''

    for bam, mRNA in zip(bamfiles, mRNAbamfiles):
        path, file = os.path.split(bam)
        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))

        dbed, bednames, bedids = read_bed(srcbedfile)
        dgff3, gff3names = read_gff3(gff3file)
        dscafflengths, dgenome, names = length_of_fasta_sequences(genomefile)
        max_id = max([int(v[3][1:]) for k, v in dbed.items()])  # get max id number # c1, c2, c3, ..., c416, ..

        if libsize:  # if it is int() (not zero)
            srnaallreads = libsize
            mrnaallreads, mrnamappedreads, mrnaunmappedreads = count_all_reads_in_bam(mRNA)
        else:
            srnaallreads, srnamappedreads, srnaunmappedreads = count_all_reads_in_bam(bam)
            mrnaallreads, mrnamappedreads, mrnaunmappedreads = count_all_reads_in_bam(mRNA)

        print(f'Limiting per base coordinate plots by file {limitbed1}')
        dbedall, bednamesall, bedidsall = read_bed(limitbedall)
        dbed1, bednames1, bedids1 = read_bed(limitbed1)
        dbed2, bednames2, bedids2 = read_bed(limitbed2)
        # target windows = [scaffold, start, end]  # limit clusters by entries in limitbedall
        targetwindows = [[dbedall[k][0], dbedall[k][1], dbedall[k][2]] for k in bednamesall]

        strandwindows, pcatable, pcacolumnnames, pcarownames, pcaindex, newclusters = \
            strandedness(bam, mRNA, targetwindows, dbed, bednames, dgff3, gff3names, dgenome, max_id, srnaallreads, mrnaallreads)

        strandedfilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.tsv' % (prefix, step, winsize, cutoff)
        write_to_bed(strandwindows, strandedfilename)
        pcatablefilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.pcatable.tsv' % (prefix, step, winsize, cutoff)
        write_table(pcatable, pcatablefilename)

        # make pca and pairwise plots
        plot_pca_and_pairwise_correlation(pcatable, pcacolumnnames, pcaindex, prefix, varexp)

        #######
        ####### reads oriented by genome, windows used limited by tartetwindows
        dclusters, clusternames = per_base_strandedness(bam, targetwindows)  # refinedwindowfilename
        plot_reads_oriented_by_genome(dclusters, clusternames, srnaallreads, dbed, bednames, numplots, prefix)

        #######
        ####### orient clusters by mRNA. Flip mRNA so always pointing in forward direction, change srna accordingly
        # these will be 'polar' variables, here bedids1 limit clusters
        dclusterspolar, clusternamespolar = per_base_strandedness_orient_by_mrna_limit(bam, strandedfilename, bedids1)
        plot_reads_oriented_by_mRNA(dclusterspolar, clusternamespolar, srnaallreads, dbed, bednames, numplots, prefix)

        #######
        ####### orient clusters by gene feature ID, do not include clusters covered by lots of mRNA (~Class II only..)
        # here bedids2 limit clusters
        dclustersgene, clusternamesgene = per_base_strandedness_orient_by_gene_limit(bam, strandedfilename, bedids2)
        plot_reads_oriented_by_gene(dclustersgene, clusternamesgene, srnaallreads, dbed, bednames, numplots, prefix)

    return


def main_denovo(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step=25, winsize=100, cutoff=100, numplots=4,
                varexp=0.9, libsize=0, srcbedfilefrompaper='', update_srcs=False):
    import os
    from natsort import natsorted, ns
    import pandas as pd

    for bam, mRNA in zip(bamfiles, mRNAbamfiles):
        path, file = os.path.split(bam)
        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))

        dbed, bednames, bedids = read_bed(srcbedfile)
        dgff3, gff3names = read_gff3(gff3file)
        dscafflengths, dgenome, names = length_of_fasta_sequences(genomefile)
        max_id = max([int(v[3][1:]) for k, v in dbed.items()])  # get max id number # c1, c2, c3, ..., c416, ..

        if libsize:  # if it is int() (not zero)
            srnaallreads = libsize
            mrnaallreads, mrnamappedreads, mrnaunmappedreads = count_all_reads_in_bam(mRNA)
        else:
            srnaallreads, srnamappedreads, srnaunmappedreads = count_all_reads_in_bam(bam)
            mrnaallreads, mrnamappedreads, mrnaunmappedreads = count_all_reads_in_bam(mRNA)


        print('Denovo cluster identification: step: %d, window size: %d, cutoff: %d' % (step, winsize, cutoff))

        windows, positions = calculate_depth_with_pysam_sliding_windows(bam, dscafflengths, step, winsize, cutoff)

        sortwindows = natsorted(['\t'.join(w) for w in windows], key=lambda y: y.lower())
        write_to_bed(sortwindows, '%s.windows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff))
        collapsedsortwindows = collapse_windows(sortwindows)
        #collapsedwindowfilename = '%s.collapsedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        print(f'Starting to refine window coordinates.')
        refinedsortwindows = refine_window_coordinates(collapsedsortwindows, bam)
        refinedwindowfilename = '%s.refinedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        write_to_bed(refinedsortwindows, refinedwindowfilename)
        targetwindows = [[win.split('\t')[0], win.split('\t')[1], win.split('\t')[2]] for win in refinedsortwindows]

        strandwindows, pcatable, pcacolumnnames, pcarownames, pcaindex, newclusters = \
            strandedness(bam, mRNA, targetwindows, dbed, bednames, dgff3, gff3names, dgenome, max_id, srnaallreads, mrnaallreads)
        ### add the new clusters to the SRC file so they can be considered in the next run.
        if srcbedfilefrompaper:
            count_new_srcs_compared_to_paper(srcbedfilefrompaper, strandwindows)
        if update_srcs:
            update_src_bed_file_with_new_srcs(srcbedfile, newclusters)
        strandedfilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.tsv' % (prefix, step, winsize, cutoff)
        write_to_bed(strandwindows, strandedfilename)
        pcatablefilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.pcatable.tsv' % (prefix, step, winsize, cutoff)
        write_table(pcatable, pcatablefilename)

        # make pca and pairwise plots
        plot_pca_and_pairwise_correlation(pcatable, pcacolumnnames, pcaindex, prefix, varexp)

        #######
        ####### reads oriented by genome
        dclusters, clusternames = per_base_strandedness(bam, targetwindows)  # refinedwindowfilename
        plot_reads_oriented_by_genome(dclusters, clusternames, srnaallreads, dbed, bednames, numplots, prefix)

        #######
        ####### orient clusters by mRNA. Flip mRNA so always pointing in forward direction, change srna accordingly
        # these will be 'polar' variables
        dclusterspolar, clusternamespolar = per_base_strandedness_orient_by_mrna(bam, strandedfilename)
        plot_reads_oriented_by_mRNA(dclusterspolar, clusternamespolar, srnaallreads, dbed, bednames, numplots, prefix)

        #######
        ####### orient clusters by gene feature ID, do not include clusters covered by lots of mRNA (~Class II only..)
        dclustersgene, clusternamesgene = per_base_strandedness_orient_by_gene(bam, strandedfilename)
        plot_reads_oriented_by_gene(dclustersgene, clusternamesgene, srnaallreads, dbed, bednames, numplots, prefix)

    return


def main(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step=25, winsize=100, cutoff=100, numplots=4,
         varexp=0.9, libsize=0, srcbedfilefrompaper='', update_srcs=False, limitbedall='', limitbed1='', limitbed2=''):
    # input is sam file # output does not include positions of 0 coverage
    # key is scafffold_position # value is depth

    if limitbedall and limitbed1 and limitbed2:
        main_limit_clusters(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step, winsize, cutoff, numplots,
                            varexp, libsize, limitbedall, limitbed1, limitbed2)
    else:
        main_denovo(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step, winsize, cutoff, numplots, varexp,
                    libsize, srcbedfilefrompaper, update_srcs)

    print('###FIN###')
