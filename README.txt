
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

