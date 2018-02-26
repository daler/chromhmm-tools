#!/usr/bin/env python

"""
`chromhmm_signalize.py`
-----------------------
Starting from bigWig or bedGraph files, this program helps construct input
files for ChromHMM's BinarizeSignal program.

ChromHMM's BinarizeSignal expects input data of this form::

    {celltype}  {chromosome}
    {mark1}     {mark2}     {mark3}
    0           4           0
    1           3           0
    2           1           9


Values are read counts in each consecutive bin of a chromosome, for each mark
(e.g., histone modifications).

If you have mapped reads, then ChromHMM's BinarizeBed is best for generating
this kind of file. However, lots of publicly availiable data are available as
bigWig files, which do not have regularly spaced bins.


The "rebin" program interpolates bigWig or bedGraph signal into regularly
spaced bins, and separates this interpolated signal by chromosome.

The "combine" program merges interpolated files into a single file suitable for
use by ChromHMM's BinarizeSignal program (and the rest of the ChromHMM
workflow).

REBIN
-----

    The config file is a tab-delimited file with 3 required and 1 optional
    field::

        {celltype} {mark} {path} [{control}]

    Each file specified by `path` will be interpolated into consecutive bins of
    size `binsize` and will be split into mulitple files based on the
    chromosomes present in the original file.  These new files will be named
    according to the scheme::

        {binned_dir}/{celltype}_{mark}_{binsize}_{chrom}.binned

    In the config file, `path` can refer to either a bigWig file or a bedGraph
    file (format is detected automatically).  If bigWig, then it will be
    converted into a bedGraph with UCSC's bigWigToBedGraph program, which must
    be on the path.  This converted file will be saved with a '.bedgraph'
    extension next to the original bigWig.  If this file already exists, it
    will be used without doing the conversion again.

COMBINE
-------

    After running `rebin`, run `combine` to combine all files of the same
    celltype and chromosome into a single "_signal" file, which can be directly
    used with ChromHMM's BinarizeSignal program.  Best practice is to re-use
    the same config file used for `rebin`.


GENERAL WORKFLOW
----------------

    * Write a config file
        (see above for format; here called `config.txt`)
    * `chromhmm-signalize.py config.txt hg19`
        (see `binned` dir for output by default)
    * `chromhmm-signalize.py config.txt binned`
        (see `signal` dir for output by default)
    * `java -jar -xm4000M ChromHMM.jar BinarizeSignal signal binarized`
        (see `binarized` for output)
    * `java -jar -xm4000M ChromHMM.jar LearnModel binarized output 9 hg19`
        (see `output` for 9-state model results)
"""

import argh
from argh import arg
import os
import logging
import pybedtools
import glob
import numpy as np
import itertools

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(msg)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def is_bigwig(fn):
    """
    checks magic number to see if we're working with a bigWig file
    """
    fh = open(fn, 'rb')
    magic = fh.read(4)
    if (magic == '&\xfc\x8f\x88') or (magic == '\x88\x8f\xfc&'):
        return True


class Config(object):
    """
    Lightweight object to parse and store config info
    """
    def __init__(self, filename):
        self.filename = filename
        self.config = []
        for line in open(filename):
            if line.startswith('#'):
                continue
            toks = line.strip().split('\t')
            if len(toks) == 1:
                continue
            celltype, mark, path = toks[:3]
            if len(toks) == 4:
                control = toks[-1]
            else:
                control = None
            self.config.append(
                (celltype, mark, path, control)
            )


@arg('configfile', help='input config file')
@arg('genome', help='assembly name')
@arg('--binsize', help='Bin size for _signal file, default is "%(default)s"')
@arg('--binned_dir', help='Output dir, default "%(default)s"')
@arg('--quiet', help='disable logging')
def rebin(configfile, genome, binsize=200, quiet=False, binned_dir='binned'):
    """
    Split bigwig/bedGraph files by chromosome, and interpolate signal into
    `binsize` bins.
    """
    if quiet:
        logger.disabled = True

    config = Config(configfile)

    genome = pybedtools.chromsizes(genome)

    if not os.path.exists(binned_dir):
        os.makedirs(binned_dir)

    for celltype, mark, path, control in config.config:
        chrom = None
        fout = None
        output_pattern = (
            '{binned_dir}/{celltype}-{mark}-{binsize}-'
            '{{chrom}}.binned'.format(**locals()))
        logger.info('{path} -> {output_pattern}'.format(**locals()))

        # convert to bedGraph if needed
        if is_bigwig(path):
            bg = path + '.bedgraph'
            if not os.path.exists(bg):
                logger.info('converting to bedgraph')
                os.system('bigWigToBedGraph %s %s' % (path, bg))
            else:
                logger.info('%s already exists, using it' % bg)
        else:
            bg = path

        bt = pybedtools.BedTool(bg)
        x = []
        y = []

        def write_interpolated_results(x, y, chrom):
            """
            interpolation and file-creation happens here
            """
            logger.info(chrom)
            filename = output_pattern.format(chrom=chrom)
            max_pos = genome[chrom][-1]
            x = np.array(x)
            y = np.array(y)
            xi = np.arange(0, max_pos, binsize)
            yi = np.interp(xi, x, y, left=-1, right=-1)

            fout = open(filename, 'w')
            fout.write('%s\t%s\n' % (celltype, chrom))
            fout.write('%s\n' % mark)
            for xii, yii in itertools.izip(xi, yi):
                fout.write('%s\n' % yii)
            fout.close()

            # try to save a little memory
            del x, y, xi, yi

        for i in bt:
            if (i.chrom != chrom) and (chrom is not None):
                write_interpolated_results(x, y, chrom)
                x = []
                y = []

            # use the midpoint of each bedgraph feature
            x.append(i.start + (i.stop - i.start) / 2)
            y.append(float(i[-1]))
            chrom = i.chrom

        # last one
        write_interpolated_results(x, y, chrom)

        if quiet:
            logger.disabled = False


@arg('configfile', help='input config file')
@arg('--binned_dir', help='Dir containing binned data, default "%(default)s"')
@arg('--binsize', help='Bin size for _signal file, default is "%(default)s"')
@arg('--combined_dir',
     help='Output dir for combined files, default is "%(default)s"')
def combine(configfile, binned_dir='binned', binsize=200,
            combined_dir='signal'):
    """
    Combine interpolated files from the "rebin" command into "_signal" files
    ready for ChromHMM's BinarizeSignal.
    """

    config = Config(configfile)
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)

    # accumluate all possible files based on the config
    files = []
    for celltype, mark, path, control in config.config:
        pattern = os.path.join(
            binned_dir,
            '{celltype}-{mark}-{binsize}-*.binned'.format(**locals()))
        files.extend(glob.glob(pattern))

    # One file for each (celltype, chrom) combo which aggregates all marks.
    def sortkey(x):
        (celltype, mark, binsize, chrom) = \
            os.path.basename(x).replace('.binned', '').split('-')
        return (celltype, chrom)

    files.sort(key=sortkey)
    for item, group in itertools.groupby(files, sortkey):
        # this group should now have all the same chromosome and celltype
        group = sorted(list(group))
        celltype, chrom = item

        logger.info('%s, %s, combining:' % (celltype, chrom))

        for g in group:
            logger.info('  ' + g)

        open_files = [open(f) for f in group]
        header_1 = [f.readline() for f in open_files]
        cells, chroms = zip(*[i.strip().split('\t') for i in header_1])
        if len(set(cells)) != 1 or list(set(cells))[0] != celltype:
            raise ValueError("expecting a single celltype, %s found" % cells)
        if len(set(chroms)) != 1 or list(set(chroms))[0] != chrom:
            raise ValueError(
                "filename chrom (%s) does not match header (%s)"
                % (chrom, chroms))

        header_2 = [f.readline() for f in open_files]
        marks = [i.strip() for i in header_2]

        filename = os.path.join(
            combined_dir,
            '{celltype}_{chrom}_{binsize}_signal.txt'.format(**locals()))
        logger.info('writing %s' % filename)
        final = open(filename, 'w')
        final.write('%s\t%s\n' % (celltype, chrom))
        final.write('\t'.join(marks) + '\n')
        for values in itertools.izip(*open_files):
            values = [v.strip() for v in values]
            final.write('\t'.join(values) + '\n')
        final.close()

if __name__ == "__main__":
    parser = argh.ArghParser()
    #parser.usage = __doc__
    parser.add_commands([rebin, combine])
    parser.dispatch()
