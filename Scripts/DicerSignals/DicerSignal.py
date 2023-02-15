######################################################################
##########   Assumes .sam files as input (w/ extension)     ##########
######################################################################



def plot_scatter(x, y, outpath, title_analysis, exitplotcount, countplot, legenlocation='upper right', prefix='NaN', xlab='Distance From Read End', ylab='Depth', figlength=10, figheight=5):
    ''' x is list of x values, y is list of y values, path is full path to output file '''
    import matplotlib.pyplot as plt
    # import math
    # import os
    # font = {'family': 'normal',
    #         'weight': 'normal',
    #         'size': 20}
    if countplot == 1:
        plt.rcParams.update({'font.size': 16})
        plt.figure(figsize=(figlength, figheight))

    # y = [math.log(i, 2) for i in y]
    plt.plot(x, y, label=prefix)
    plt.grid(True)
    plt.title(title_analysis)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend(loc=legenlocation, fontsize=10)
    plt.savefig(outpath)
    if countplot == exitplotcount:  # 2 = len(infiles), infiles passed to phasing function
        plt.show()
        plt.close()


def ping_pong(infiles_f, infiles_r, lower=0, upper=20, analysis='PingPong55'):
    ''' Assumes first forward and reverse files are from the same sam/bam file '''
    ''' Main difference from phasing analysis is comparing forward to reverse reads '''
    import os
    from natsort import natsorted, ns

    end1 = 2  # access 5' end information
    end2 = 2  # access 5' end information
    print('lower: %d, upper: %d, orientation not applicable' % (lower, upper))

    allsignals = []
    alldicer = []
    for forward, reverse in zip(infiles_f, infiles_r):
        print(forward)
        print(reverse)
        df, dr = {}, {}
        namesf, namesr = [], []
        with open(forward, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                df['_'.join(s[0:2])] = s
                namesf.append('_'.join(s[0:2]))
        with open(reverse, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                dr['_'.join(s[0:2])] = s
                namesr.append('_'.join(s[0:2]))

        dicer_signal_location = {}  # record locations and strength of dicer signal at position upper - 3 (should be - 2 if count was 1 based)
        # only need to compare forward to reverse, not reverse to forward. This would be redundant
        signals = [0] * abs(upper - lower)
        count = 0  # keep track of element in list signals that i should add too
        for x in range(lower, upper):  # range excludes the last value of 50 so use 51 # distance from 3' end
            minimums = []
            for name in namesf:  # basically for each name in forward reads
                s = df[name]  # list of elements of line previously called s
                if int(s[end1]) != 0:  # for every position with a 5' end in forward read orientation
                    # now try to compare to 5' ends in reverse file
                    Mi = int(s[end1])  # get 5' counts from current position for ping pong
                    try:
                        dr['_'.join([s[0], str(int(s[1]) + x)])]  # is same scaffold_position is present in reverse?
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in reverse file
                    else:  # get 5' counts from i+x position
                        Nix = int(dr['_'.join([s[0], str(int(s[1]) + x)])][end2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))  # min(Mi, Nix) # Mi * Nix # assume that reads must be paired 1/1 ratio
                    if (count == upper - 3) and (min(Mi, Nix) != 0):  # count is zero based so it should be -2 normally for dicer products
                        dicer_signal_location[name] = dicer_signal_location.setdefault(name, 0) + min(Mi, Nix)
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)
        alldicer.append(dicer_signal_location)

    outfilename = 'PingPong.FvsR.55'
    # add name of files to output
    output = []
    countplot = 1
    for fullpath, signal in zip(infiles_f, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'pdf'])), analysis, len(allsignals), countplot, 'upper left', '_'.join(file.split('.')[:4]))
        countplot += 1

    # output data
    path, file = os.path.split(infiles_f[0])
    outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'out']))
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

    for fullpath, dicer_dictionary in zip(infiles_f, alldicer):
        path, file = os.path.split(fullpath)
        # sort dictionary by value and output. output is position info as key and sorted by value # of 5-5
        # {k: v for k, v in sorted(x.items(), key=lambda item: item[1])}
        output = [f'{k}\t{v}' for k, v in sorted(dicer_dictionary.items(), key=lambda item: item[1], reverse=True)]
        outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'DicerSignal.SortStrength.out']))
        print(f'Writing to output file: {outfile}')
        with open(outfile, 'w') as OUT:
            OUT.write('\n'.join(output))
        output = [f'{k}\t{v}' for k, v in natsorted(dicer_dictionary.items(), key=lambda item: item[0])]
        outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'DicerSignal.SortPosition.out']))
        print(f'Writing to output file: {outfile}')
        with open(outfile, 'w') as OUT:
            OUT.write('\n'.join(output))


def get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate):
    import re
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate - 1  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'N', 'D', 'X', '=']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return rightmost_position


def count53(filelist, maxdepth, orientation='forward'):
    from natsort import natsorted, ns

    outnames = []
    for f in filelist:
        print('Working on: %s' % (f))
        # filelist = ['D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.R.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.R.sort.sam']
        if maxdepth > 0:  # if maxdepth specified chance output file name
            outfile = f[:-len('.sam')] + '.53.MaxDepth.tsv'
        else:
            outfile = f[:-len('.sam')] + '.53.tsv'
        outnames.append(outfile)

        d5 = {}
        d3 = {}
        with open(f, 'r') as FILE:
            for line in FILE:
                line = line.strip()
                if line[0] == '@':
                    pass
                else:
                    leftmost_coordinate = int(line.split('\t')[3])
                    CIGAR = line.strip().split('\t')[5]
                    if orientation == 'forward':
                        scaff_pos5 = '_'.join(line.split('\t')[2:4])  # scaffold name + 5' position
                        # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                        # scaff_pos3 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 3' position
                        scaff_pos3 = '_'.join([line.split('\t')[2]] + [str(get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate))])  # scaffold name + 3' position
                    elif orientation == 'reverse':  # switch 5' and 3' values for reverse reads
                        scaff_pos3 = '_'.join(line.split('\t')[2:4])  # scaffold name + 3' position
                        # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                        # scaff_pos5 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 5' position
                        scaff_pos5 = '_'.join([line.split('\t')[2]] + [str(get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate))])
                    # could have done below two lines, but doesn't keep order
                    # d5[scaff_pos5].get(scaff_pos5, 0) + 1  # kinda slow # counts number of 5' positions for all scaffolds
                    # d3[scaff_pos3].get(scaff_pos3, 0) + 1  # kinda slow # counts number of 3' positions for all scaffolds
                    # count 5' ends at position
                    try:
                        d5[scaff_pos5]
                    except:
                        d5[scaff_pos5] = 1  # initiate with 1 because we want to count the first
                    else:
                        d5[scaff_pos5] += 1
                    try:  # to maintain the same keys between d5 and d3
                        d5[scaff_pos3]
                    except:  # if there is no scaff_pos3 entry make it zero
                        d5[scaff_pos3] = 0

                    # count 3' ends at each position
                    try:
                        d3[scaff_pos3]
                    except:
                        d3[scaff_pos3] = 1  # initiate  with 1 because we want to count the first
                    else:
                        d3[scaff_pos3] += 1
                    try:  # to maintain the same keys between d5 and d3
                        d3[scaff_pos5]
                    except:  # if there is no scaff_pos5 entry in d3 make it zero
                        d3[scaff_pos5] = 0
        # sort key values naturally # scaffolds may not be in same order as sam but coordinates within each scaffold should be ordered
        names = list(d5.keys())  # assumes d5 and d3 have same key values
        sortnames = natsorted(names, key=lambda y: y.lower())  # or # natsorted(x, alg=ns.IGNORECASE)  # or alg=ns.IC
        # print('First sorted names:')
        # print(sortnames[:4])

        with open(outfile, 'w') as OUT:
            output = []
            for name in sortnames:
                scaffold, position = '_'.join(name.split('_')[:-1]), name.split('_')[-1]
                if maxdepth == 0:  # then do not control depth # proceed without changing depth values
                    output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
                else:  # enforce specified maxdepth value
                    if d5[name] > maxdepth:
                        d5[name] = maxdepth
                    if d3[name] > maxdepth:
                        d3[name] = maxdepth
                    output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
            OUT.write('\n'.join(output))
        print('Output file of 5\' and 3\' counts: %s\n' % (outfile))
    return(outnames)


def main(samfilelist_f, samfilelist_r, seqlength=23, maxdepth=0):  # 0 maxdepth will not limit depth > 0 will limit to max depth
    print('##################################################')
    print('Count forward and reverse reads')
    infiles_f = count53(samfilelist_f, maxdepth, orientation='forward')  # infiles = scaffold\tposition\t#of5'ends\t#of3'ends
    infiles_r = count53(samfilelist_r, maxdepth, orientation='reverse')
    print('##################################################')
    print('Starting ping pong analysis')
    ping_pong(infiles_f, infiles_r, 0, seqlength)  # compare 5' forward to 5' reverse reads and vice versa
    print('Done with ping pong analysis')
    print('##################################################')

