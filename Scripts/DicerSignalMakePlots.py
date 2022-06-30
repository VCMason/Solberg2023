

def plot_scatter_and_barplot_subplots(pingpongplotdata, barplotdata, rowtitles, outfilepath, maxy=False):
    import matplotlib.pyplot as plt
    import math

    numsubplots = len(pingpongplotdata)
    fig, axs = plt.subplots(numsubplots, 2, figsize=(10, 10), gridspec_kw={'width_ratios': [1, 2]})  # sharex=False
    axs[0, 0].set_title('Ping-Pong Signal')
    axs[0, 1].set_title('Dicer Signal Strength')

    if maxy:
        allpingy, allbary = [], []
        for i, t in enumerate(zip(pingpongplotdata, barplotdata, rowtitles)):
            scatterdata, bardata, rowname = t
            # ax[row(up/down), column(left/right)]
            pingx, pingy, pingfile = scatterdata
            barx, bary, barfile = bardata
            # rowname = '_'.join(barfile.split('_')[:2])
            pingy = [float(y) for y in pingy]
            bary = [float(y) for y in bary]
            allpingy = allpingy + pingy
            allbary = allbary + bary
        allpingymax, allbarymax = max(allpingy), max(allbary)
        # logallpingymax, logallbarymax = round(math.log(allpingymax+1, 2), 3), round(math.log(allbarymax+1, 2), 3)

    for i, t in enumerate(zip(pingpongplotdata, barplotdata, rowtitles)):
        scatterdata, bardata, rowname = t
        # ax[row(up/down), column(left/right)]
        pingx, pingy, pingfile = scatterdata
        barx, bary, barfile = bardata
        # rowname = '_'.join(barfile.split('_')[:2])
        pingx = [int(x) for x in pingx]
        pingy = [float(y) for y in pingy]
        bary = [float(y) for y in bary]
        pingymin, pingymax = round(min(pingy), 3), round(max(pingy), 3)
        barymin, barymax = round(min(bary), 3), round(max(bary), 3)
        axs[i, 0].plot(pingx, pingy)  # color=color1, label=label1 # pingx is list of ints, pingy is list of ints
        axs[i, 1].bar(barx, bary)  # barx is list of strings, bary is list of integers
        if maxy:
            # logpingy = [round(math.log(y+1, 2), 3) for y in pingy]
            # logbary = [round(math.log(y+1, 2), 3) for y in bary]
            # logpingymax, logbarymax = round(math.log(pingymax+1, 2), 3), round(math.log(barymax+1, 2), 3)
            # logpingymin, logbarymin = round(math.log(pingymin+1, 2), 3), round(math.log(barymin+1, 2), 3)
            axs[i, 0].set_ylim([0, allpingymax])  # pingymin
            axs[i, 1].set_ylim([0, allbarymax])  # barymin
            axs[i, 0].set_yticks([pingymax])
            axs[i, 1].set_yticks([barymax])
            # axs[i, 0].set_yticks([logpingymin, logpingymax, logallpingymax])
            # axs[i, 1].set_yticks([logbarymin, logbarymax, logallbarymax])
            # if allpingymax in pingy:
            #     axs[i, 0].set_yticks([logpingymin, logallpingymax])
            # else:
            #     axs[i, 0].set_yticks([logpingymax, logallpingymax])
            # if allbarymax in bary:
            #     axs[i, 1].set_yticks([logbarymin, logallbarymax])
            # else:
            #     axs[i, 1].set_yticks([logbarymax, logallbarymax])
        else:
            axs[i, 0].set_ylim([pingymin, pingymax])
            axs[i, 1].set_ylim([pingymin, barymax])
            axs[i, 0].set_yticks([pingymin, pingymax])
            axs[i, 1].set_yticks([barymin, barymax])
        axs[i, 0].set_ylabel(rowname, rotation=90)
        # axs[i, 0].yaxis.set_label_position('left')
        axs[i, 1].tick_params('x', labelrotation=45, labelsize=8)
        if i+1 != len(rowtitles):
            axs[i, 0].get_xaxis().set_visible(False)
            axs[i, 1].get_xaxis().set_visible(False)

    # for ax in axs.flat:
    #     ax.set(xlabel='x-label', ylabel='y-label')
    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()

    plt.savefig(outfilepath)
    plt.show()
    plt.close()

    return


def barplot(clusters, strengths, outpath, title='', xlab='Ptiwi08 Bound sRNA Cluster', ylab='Dicer Signal Strength', figlength=10, figheight=5):
    # strengths is list of integers representing height of bar plot
    # clusters is list of cluster names (strings), representing bar plot categories
    import matplotlib.pyplot as plt
    import os

    p, f = os.path.split(outpath)
    title = f[:-len('.F.sort.53.PingPong.FvsR.55.DicerSignal.SortPosition.pdf')]
    #general figure setup
    plt.rcParams.update({'font.size': 12})
    plt.rcParams.update({'figure.autolayout': True})  # automatically adjust fig so x axis labels are visible
    plt.figure(figsize=(figlength, figheight))
    plt.xticks(rotation=45)  # rotation='vertical' # , fontsize=12 # rotation=45, ha='right'
    # plot data and axis labels
    plt.bar(clusters, strengths)
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(outpath)
    plt.show()
    plt.close()
    return


def summarize_dicer_signal_per_cluster(bed, dpp, totalcount):

    dclusterstrength = {}
    names = []
    for i in bed:
        scaffold, start, end = i[0], int(i[1]), int(i[2])
        key = i[3]  # cluster id
        names.append(key)
        dclusterstrength[key] = 0
        for position, strength in dpp.items():
            ppscaff = '_'.join(position.split('_')[:-1])
            ppcoord = int(position.split('_')[-1])
            if (ppscaff == scaffold) and (ppcoord >= start) and (ppcoord <= end):
                dclusterstrength[key] = dclusterstrength.setdefault(key, 0) + int(strength)
    # normalize by total reads mapping inside bed features
    dclusterstrengthnorm = {key: round((val / totalcount)*100, 5) for key, val in dclusterstrength.items()}

    return dclusterstrengthnorm, names


def collect_dicer_strength(pingpingfile):
    # scaffold51_8_23317  2
    # scaffold51_8_23341  2
    # scaffold51_8_23348  1
    # scaffold51_8_23351  6
    # scaffold51_8_23353  19
    with open(pingpingfile, 'r') as FILE:
        dpingpong = {line.strip().split('\t')[0]: int(line.strip().split('\t')[1]) for line in FILE if line.strip() != ''}

    return dpingpong


def read_ping_pong_file(pingpongfile, totalcounts):
    # import os
    # path, f = os.path.split(pingpongfile)

    pingpongplotdata = []
    lower, upper = 0+1, 23+1  # add 1 so x axis is 1 based in future plot
    with open(pingpongfile, 'r') as FILE:
        for line, totalcount in zip(FILE, totalcounts):
            filename = line.strip().split('\t')[0]
            signal = line.strip().split('\t')[1:]
            x = list(range(lower, upper))
            y = [round((int(sig) / totalcount)*100, 4) for sig in signal]
            # title = '_'.join(filename.split('.')[:4])
            # outpath = os.path.join(path, '.'.join(filename.split('.')[:-1] + ['pdf']))
            pingpongplotdata.append([x, y, filename])

    return pingpongplotdata


def total_read_count_from_all_windows(bamfilename, bed):
    """bed is list of lists of each line and tab delimited element in bed file, header excluded"""
    """bamfilename is full path to bam file"""
    import pysam

    totalcount = 0
    samfile = pysam.AlignmentFile(bamfilename, "rb")
    for line in bed:
        scaff, start, end = line[0], int(line[1]), int(line[2])
        for read in samfile.fetch(scaff, start, end):
            totalcount += 1

    return totalcount


def read_bed(bedfile):
    # scaffold51_8	23301	23738   cluster8
    with open(bedfile, 'r') as FILE:
        # scaffold = line.strip().split('\t')[0]
        # start = int(line.strip().split('\t')[1])
        # end = int(line.strip().split('\t')[2])
        # id = line.strip().split('\t')[3]
        bed = [line.strip().split('\t') for line in FILE if line.strip() != '']

    return bed


def main(pingpongfile, dicersignalfiles, bamfiles, rowtitles, bedfile, maxy=False):
    import os
    # pingpongfile is full path to file containing all summarized ping pong results, one line per file/sample
    # dicersignalfiles is list of files, list of files should match the order in ping pong file!
    # bedfile has at least four tab delimited columns. Scaffold start   end id
    bed = read_bed(bedfile)

    totalcounts = []  # list of integers, total counts within bed features per file
    for bamfilename in bamfiles:
        totalcount = total_read_count_from_all_windows(bamfilename, bed)
        totalcounts.append(totalcount)
        print(f'Total reads inside bed features: {totalcount} for file: {bamfilename}')

    pingpongplotdata = read_ping_pong_file(pingpongfile, totalcounts)

    barplotdata = []
    for filepath, totalcount in zip(dicersignalfiles, totalcounts):
        dpingpong = collect_dicer_strength(filepath)
        dclusterstrength, names = summarize_dicer_signal_per_cluster(bed, dpingpong, totalcount)

        path, filename = os.path.split(filepath)
        outpath = os.path.join(path, '.'.join(filename.split('.')[:-1] + ['pdf']))
        strengths = [dclusterstrength[n] for n in names]
        barplot(names, strengths, outpath)
        barplotdata.append([names, strengths, filename])

    path, filename = os.path.split(pingpongfile)
    if maxy:
        outfile = os.path.join(path, 'PingPong.DicerStrength.plots.MaxYAxisFixed.pdf')
    else:
        outfile = os.path.join(path, 'PingPong.DicerStrength.plots.pdf')
    plot_scatter_and_barplot_subplots(pingpongplotdata, barplotdata, rowtitles, outfile, maxy)
