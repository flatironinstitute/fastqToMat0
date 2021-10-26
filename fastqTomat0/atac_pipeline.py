import pandas as pd
import pandas.errors as pde

from pandas.errors import EmptyDataError
import subprocess
import os
import sys
import io
import logging

CPU_PER_TASK = 4

BWM_OUT_FILE = "aligned.sam"
BAM_UNSORTED_OUT_FILE = "aligned.bam"
BAM_FIXMATES_OUT_FILE = "aligned.fixmates.bam"
BAM_TEMP_OUT_FILE = "aligned.temp.bam"
BAM_SORTED_OUT_FILE = "aligned.sorted.bam"
BAM_DEDUPED_OUT_FILE = "aligned.sorted.dedupe.bam"

BWM_CMD = ["bwa-mem2", "mem", "-t"]
MACS_CMD = ["macs3", "callpeak", "-f", "BAMPE"]

SAMTOOLS_SORT_CMD = ["samtools", "sort", "-@"]
SAMTOOLS_FIXMATES_CMD = ["samtools", "fixmate", "-m", "-r", "-@"]
SAMTOOLS_MARKDUP_CMD = ["samtools", "markdup", "-r", "-@"]
SAMTOOLS_VIEW_CMD = ["samtools", "view", "-@"]
SAMTOOLS_INDEX_CMD = ["samtools", "index", "-@"]

CHROM_COL = 'chrom'
MULTIINTER_FULL_EXT = "_multiinter.tsv"
MULTIINTER_BED_EXT = ".bed"
MULTIINTER_COLNAMES = [CHROM_COL, 'start', 'end', 'num', 'list']
MULTIINTER_BED_COLS = [CHROM_COL, 'start', 'end']

BEDTOOLS_MULTIINTER = ['bedtools', 'multiinter', '-i']
BEDTOOLS_MERGE = ['bedtools', 'merge', '-i']
BEDTOOLS_SORT = ['bedtools', 'sort', '-i']

COUNTS_PER_CHR_CMD = ["samtools", "idxstats"]
VIEW_MAPPED_CMD = ["samtools", "view", "-f66"]
FRAGMENT_LEN_CMD = ["cut", "-f", "9"]

MACS_FILES = ["_peaks.narrowPeak", "_control_lambda.bdg", "_peaks.xls", "_summits.bed", "_treat_pileup.bdg"]

CHR_COL = "Chromosome"
COUNT_COL = "Counts"
FLEN_COL = "Fragment_Length"

COUNT_FILE_EXT = "_count.tsv"
COUNT_COLUMNS = [CHROM_COL, 'start', 'end', 'count', 'num_bases_covered', 'length', 'fraction_covered']
BEDTOOLS_COUNT = ['bedtools', 'coverage', '-sorted']

logger = logging.Logger('ATAC')
logger_handler = logging.StreamHandler(sys.stderr)
logger_handler.setFormatter(logging.Formatter('%(asctime)-15s %(message)s %(sname)s: %(cmd)s'))
logger.addHandler(logger_handler)

def process_atac_to_bed(bwa_index, fastq_r1, fastq_r2, genome_size, out_path=".", sample_name=None,
                        n_threads=CPU_PER_TASK, debug=False, nomodel=False, extsize=None):
    """
    Process FASTQ files into an aligned, deduplicated BAM file and a set of MACS peaks

    :param bwa_index: Path to BWA-MEM2 index
    :type bwa_index: str
    :param fastq_r1: Path to FASTQ read1 file(s).
    :type fastq_r1: str, list(str)
    :param fastq_r2: Path to FASTQ read2 file(s).
    :type fastq_r2: str, list(str)
    :param genome_size: Genome size in bases. Pass as a string (e.g. "2.7e9" for hs, "1.87e9" for mm, etc)
    :type genome_size: str
    :param out_path: Location for output (and temp) files, defaults to "."
    :type out_path: str, optional
    :param sample_name: Sample name, defaults to None
    :type sample_name: str, optional
    :param n_threads: Number of threads for called programs, defaults to 4
    :type n_threads: int, optional

    :return: Returns file names for bam & macs files and a tuple of summary stats
    :rtype: str, str, (int, int, pd.DataFrame, pd.DataFrame)
    """

    _fql_1, _fql_2 = isinstance(fastq_r1, (tuple, list)), isinstance(fastq_r2, (tuple, list))

    if out_path is "." and sample_name is not None:
        out_path = os.path.join(os.path.abspath(os.path.expanduser(out_path)), sample_name)
    else:
        out_path = os.path.abspath(os.path.expanduser(out_path))
    
    os.makedirs(out_path, exist_ok=True)

    if debug:
        logger.setLevel('DEBUG')

    if _fql_1 != _fql_2:
        raise ValueError("Pass both R1 and R2 as a list or as a string")
    if _fql_1 and _fql_2:
        _multi_fqs = True
    elif (_fql_1 and _fql_2) and (len(fastq_r1) != len(fastq_r2)):
        raise ValueError("Pass the same number of R1 and R2 fastq file numbers")
    else:
        _multi_fqs = False

    _created_files = []

    try:
        if _multi_fqs:
            bwa_aligned = [_align_atac_experiment(bwa_index, r1, r2, out_path=out_path, sample_name=sample_name,
                                                  n_threads=n_threads, append_to_existing=i > 0)
                           for i, (r1, r2) in enumerate(zip(fastq_r1, fastq_r2))][0]
        else:
            bwa_aligned = _align_atac_experiment(bwa_index, fastq_r1, fastq_r2, out_path=out_path,
                                                 sample_name=sample_name, n_threads=n_threads)

        _created_files.append(bwa_aligned)
        bam_file = _call_sam_to_bam(bwa_aligned, out_path, sample=sample_name, sort=True, sort_by_name=True,
                                    n_threads=n_threads)

        _remove_file(bwa_aligned)
        _created_files.append(bam_file)

        n_aligned_w_dups = _call_count_aligned(bam_file, sample=sample_name)
        bam_deduped = _samtools_remove_dups(bam_file, out_path=out_path, n_threads=n_threads, sample=sample_name)
        _remove_file(bam_file)

        n_aligned_without_dups = _call_count_aligned(bam_deduped, sample=sample_name)
        macs_file = _call_macs3(bam_deduped, out_path, genome_size, sample=sample_name,
                                nomodel=nomodel, extsize=extsize)

    finally:
        for f in _created_files:
            _remove_file(f)

    return sample_name, bam_deduped, macs_file, (n_aligned_without_dups, n_aligned_w_dups,
                                                 _get_chr_alignment_counts(bam_deduped, sample=sample_name),
                                                 _get_fragment_length(bam_deduped))


def _align_atac_experiment(bwa_index, fastq_r1, fastq_r2, out_path=None, skip_if_exists=True, n_threads=CPU_PER_TASK,
                           append_to_existing=False, sample_name=None):

    out_path = "." if out_path is None else out_path
    out_file = os.path.join(out_path, BWM_OUT_FILE)

    if skip_if_exists and os.path.exists(out_file) and append_to_existing is False:
        return out_file

    if sample_name is not None:
        print("Aligning {s} ({f})".format(s=sample_name, f=bwa_index))

    b_cmd = BWM_CMD + [str(n_threads), bwa_index, fastq_r1, fastq_r2]

    out_mode = "w+" if append_to_existing else "w"
    with open(out_file, mode=out_mode) as out_fh:
        logger.debug("[BWM ALIGN]", extra={'sname': sample_name, 'cmd': " ".join(b_cmd)})
        proc = subprocess.run(b_cmd, stdout=out_fh, stderr=subprocess.DEVNULL)

    if proc.returncode != 0:
        _remove_file(out_file)

        _msg = "bwa mem returned {c} [{cmd}]".format(
            c=proc.returncode, cmd=" ".join(b_cmd))
        raise RuntimeError(_msg)

    return out_file


def _call_sam_to_bam(file_name, out_path, sample=None, sort=True, sort_by_name=False, n_threads=CPU_PER_TASK):

    out_file = os.path.join(out_path, BAM_UNSORTED_OUT_FILE if not sort else BAM_SORTED_OUT_FILE)
    print("Converting {f} to BAM".format(f=file_name))

    if sort:
        return _call_sort_bam(file_name, out_file=out_file, by_name=sort_by_name, n_threads=n_threads, sample=sample)
    else:
        sb_cmd = SAMTOOLS_VIEW_CMD + ["-S", "-b", file_name]

        with open(out_file, mode="w") as out_fh:
            logger.debug("[SAMTOOLS VIEW]", extra={'sname': sample, 'cmd': " ".join(sb_cmd)})
            proc = subprocess.run(sb_cmd, stdout=out_fh, stderr=subprocess.DEVNULL)

        if proc.returncode != 0:
            _remove_file(out_file)

            _msg = "{s} returned {c} [{cmd}]".format(s=sample, c=proc.returncode, cmd=" ".join(sb_cmd))
            raise RuntimeError(_msg)

        return out_file


def _call_sort_bam(file_name, out_path=None, out_file=None, sample=None, by_name=False, n_threads=CPU_PER_TASK):

    out_file = os.path.join(out_path, BAM_SORTED_OUT_FILE) if out_file is None else out_file
    print("Sorting {f}".format(f=file_name))

    sb_cmd = SAMTOOLS_SORT_CMD + [str(n_threads)]

    if by_name:
        sb_cmd = sb_cmd + ["-n"]

    sb_cmd = sb_cmd + ["-o", out_file, file_name]

    logger.debug("[SAMTOOLS SORT]", extra={'sname': sample, 'cmd': " ".join(sb_cmd)})
    proc = subprocess.run(sb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if proc.returncode != 0:
        _remove_file(out_file)

        _msg = "{s} returned {c} [{cmd}]".format(s=sample, c=proc.returncode, cmd=" ".join(sb_cmd))
        raise RuntimeError(_msg)

    # Also index the sorted BAM file just in case
    if not by_name:
        _call_bam_index(out_file, sample=sample)

    return out_file


def _call_bam_index(bam_file, sample=None, n_threads=CPU_PER_TASK):

    # Also index the sorted BAM file just in case
    sb_idx_cmd = SAMTOOLS_INDEX_CMD + [str(n_threads), bam_file]

    logger.debug("[SAMTOOLS INDEX]", extra={'sname': sample, 'cmd': " ".join(sb_idx_cmd)})
    proc = subprocess.run(sb_idx_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if proc.returncode != 0:
        _remove_file(bam_file + ".bai")

        _msg = "{s} returned {c} [{cmd}]".format(s=sample, c=proc.returncode, cmd=" ".join(sb_idx_cmd))
        raise RuntimeError(_msg)

    return bam_file

def _samtools_remove_dups(bam_file, out_path, resort_by_name=False, n_threads=CPU_PER_TASK, sample=None):

    out_file = os.path.join(out_path, BAM_DEDUPED_OUT_FILE)
    output_nsort_file_name, output_fixmates_file_name = None, None

    try:

        # Sort by name if necessary
        if resort_by_name:
            output_nsort_file_name = _call_sort_bam(bam_file,
                                                    out_file=os.path.join(out_path, BAM_TEMP_OUT_FILE), by_name=True)

        # Run fixmates to remove secondaries and add ms values
        output_fixmates_file_name = os.path.join(out_path, BAM_FIXMATES_OUT_FILE)

        print("Running samtools fixmates on {f}".format(f=output_nsort_file_name if resort_by_name else bam_file))
        samtools_fixmates_cmd = SAMTOOLS_FIXMATES_CMD + [str(n_threads)]
        samtools_fixmates_cmd = samtools_fixmates_cmd + [output_nsort_file_name if resort_by_name else bam_file]
        samtools_fixmates_cmd = samtools_fixmates_cmd + [output_fixmates_file_name]

        logger.debug("[SAMTOOLS FIXMATES]", extra={'sname': sample, 'cmd': " ".join(samtools_fixmates_cmd)})
        proc = subprocess.run(samtools_fixmates_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if proc.returncode != 0:
            _msg = "samtools fixmates returned {c} [{cmd}]".format(c=proc.returncode, cmd=" ".join(samtools_fixmates_cmd))
            raise RuntimeError(_msg)

        _remove_file(output_nsort_file_name)

        # Resort by coordinates
        output_resort_file_name = os.path.join(out_path, BAM_TEMP_OUT_FILE)
        output_resort_file_name = _call_sort_bam(output_fixmates_file_name, out_file=output_resort_file_name)
        _remove_file(output_fixmates_file_name)

        # Remove duplicates with markdup
        print("Running samtools markdup on {n}: {f}".format(n=sample, f=output_resort_file_name))
        samtools_markdup_cmd = SAMTOOLS_MARKDUP_CMD + [str(n_threads), output_resort_file_name, out_file]

        logger.debug("[SAMTOOLS MARKDUP]", extra={'sname': sample, 'cmd': " ".join(samtools_markdup_cmd)})
        proc = subprocess.run(samtools_markdup_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if proc.returncode != 0:
            _msg = "samtools markdup returned {c} [{cmd}]".format(c=proc.returncode, cmd=" ".join(samtools_markdup_cmd))
            raise RuntimeError(_msg)

        _call_bam_index(out_file)

    finally:
        _remove_file(output_resort_file_name)
        _remove_file(output_resort_file_name + ".bai")
        _remove_file(output_nsort_file_name)
        _remove_file(output_fixmates_file_name)

    return out_file


def _call_count_aligned(file_name, sample=None, n_threads=CPU_PER_TASK):

    sb_cmd = SAMTOOLS_VIEW_CMD + [str(n_threads), "-c", "-F", "260", file_name]

    print("Counting aligned reads in {f}".format(f=file_name))

    logger.debug("[SAMTOOLS COUNT]", extra={'sname': sample, 'cmd': " ".join(sb_cmd)})
    proc = subprocess.run(sb_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    if proc.returncode != 0:
        _msg = "{s} returned {c} [{cmd}]".format(s=sample, c=proc.returncode, cmd=" ".join(sb_cmd))
        print(_msg)

        return 0

    return int(proc.stdout.decode('utf-8'))


def _get_chr_alignment_counts(bam_file, sample=None):

    qc_cmd = COUNTS_PER_CHR_CMD + [bam_file]

    logger.debug("[SAMTOOLS IDXSTATS]", extra={'sname': sample, 'cmd': " ".join(qc_cmd)})
    qc_proc = subprocess.run(qc_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    if qc_proc.returncode != 0:
        _msg = "{s} returned {c} [{cmd}]".format(s=bam_file, c=qc_proc.returncode, cmd=" ".join(qc_cmd))
        raise RuntimeError(_msg)

    try:
        df = pd.read_csv(io.StringIO(qc_proc.stdout.decode("utf-8")), sep="\t", header=None)
    except EmptyDataError:
        _msg = "{s} returned {c} [{cmd}]".format(s=bam_file, c=qc_proc.returncode, cmd=" ".join(qc_cmd))
        raise RuntimeError(_msg)

    df = df.iloc[:, [0, 2]].copy()
    df.columns = [CHR_COL, COUNT_COL]
    df.set_index(CHR_COL, inplace=True)

    return df


def _get_fragment_length(bam_file, max_len=2000, binwidth=10, sample=None):

    frag_len_cmd = VIEW_MAPPED_CMD + [bam_file, "|"] + FRAGMENT_LEN_CMD

    logger.debug("[SAMTOOLS VIEW FRAG LEN]", extra={'sname': sample, 'cmd': " ".join(frag_len_cmd)})
    flen_proc = subprocess.run(" ".join(frag_len_cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    df = pd.read_csv(io.StringIO(flen_proc.stdout.decode("utf-8")), names=[FLEN_COL])
    c = pd.cut(df[FLEN_COL], pd.interval_range(start=0, end=max_len, freq=binwidth)).value_counts()
    c.index = c.index.map(lambda x: x.left).astype(int)
    c = c.reindex(range(0, max_len, binwidth))

    return c


def _remove_file(file, msg=False):

    if file is None:
        return
    else:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass
        else:
            logger.debug("[REMOVE FILE]", extra={'sname': "rm", 'cmd': file})
    if msg:
        print("Removed file {f}".format(f=file))


def _call_macs3(file_name, out_path, genome_size, sample=None, nomodel=False, extsize=None):

    sample = "NA" if sample is None else sample
    out_file = os.path.join(out_path, sample + MACS_FILES[0])

    if os.path.exists(out_file):
        return out_file

    print("Calling peaks with MACS3 from {f}".format(f=file_name))

    m3_cmd = MACS_CMD + ["-t", file_name, "-g", genome_size, "-n", sample, "-B", "-q", "0.01",
                         "--outdir", out_path]

    if nomodel:
        m3_cmd = m3_cmd + ["--nomodel"] if extsize is None else ["--nomodel", str(extsize)]

    logger.debug("[MACS3]", extra={'sname': sample, 'cmd': " ".join(m3_cmd)})
    proc = subprocess.run(m3_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if proc.returncode != 0:
        _msg = "{s} returned {c} [{cmd}]".format(s=sample, c=proc.returncode, cmd=" ".join(m3_cmd))

        for fn in MACS_FILES:
            _remove_file(os.path.join(out_path, sample + fn))

        raise RuntimeError(_msg)

    return out_file


def _bedtools_multiinter(in_files, out_path, sample="NA", bed_sort_file=None):

    out_file = os.path.join(out_path, sample + MULTIINTER_BED_EXT)

    if os.path.exists(out_file):
        return out_file

    _keep_file = False
    try:

        try:
            print("Running bedtools multiinter for {f}".format(f=sample))

            multiinter_cmd = BEDTOOLS_MULTIINTER + [i[1] for i in in_files]
            logger.debug("[BEDTOOLS MULTIINTER]", extra={'sname': sample, 'cmd': " ".join(multiinter_cmd)})
            multiinter_proc = subprocess.run(multiinter_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            multiinter_df = pd.read_csv(io.StringIO(multiinter_proc.stdout.decode("utf-8")), sep="\t", header=None)
            multiinter_df.columns = MULTIINTER_COLNAMES + [i[0] for i in in_files]

            multiinter_df.to_csv(os.path.join(out_path, sample + MULTIINTER_FULL_EXT), sep="\t")

        except (RuntimeError, pde.EmptyDataError):
            print(io.StringIO(multiinter_proc.stdout.decode("utf-8")).read(1000))
            print(io.StringIO(multiinter_proc.stderr.decode("utf-8")).getvalue())
            raise 

        try:
            print("Combining peaks into single BED file for {f}".format(f=sample))

            multiinter_df[MULTIINTER_BED_COLS].to_csv(out_file, sep="\t", header=None, index=False)

            merge_proc = subprocess.run(BEDTOOLS_MERGE + [out_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.debug("[BEDTOOLS MERGE]", extra={'sname': sample, 'cmd': " ".join(BEDTOOLS_MERGE + [out_file])})

            merge_df = pd.read_csv(io.StringIO(merge_proc.stdout.decode("utf-8")), sep="\t", header=None)
            merge_df.to_csv(out_file, sep="\t", index=False, header=None)

            if bed_sort_file is not None:

                print("Sorting output BED file")

                sort_cmd = BEDTOOLS_SORT + [out_file, "-g", bed_sort_file]
                logger.debug("[BEDTOOLS SORT]", extra={'sname': sample, 'cmd': " ".join(sort_cmd)})
                
                sort_proc = subprocess.run(sort_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.PIPE)
                sort_df = pd.read_csv(io.StringIO(sort_proc.stdout.decode("utf-8")), sep="\t", header=None)
                sort_df.to_csv(out_file, sep="\t", index=False, header=None)

            _keep_file = True

            return out_file

        except (RuntimeError, pde.EmptyDataError):
            print(" ".join(BEDTOOLS_MERGE + [out_file]))
            print(io.StringIO(merge_proc.stdout.decode("utf-8")).read(1000))
            print(io.StringIO(merge_proc.stderr.decode("utf-8")).getvalue())
            raise 

    finally:
        if not _keep_file:
            _remove_file(out_file)


def _bedtools_get_cov_count(in_files, bed_file, out_path, sample="NA"):
    """

    Creates a count table for reads that overlap a BED file intervals out of an arbitrary number of samples

    :param in_files: A list of (sample_id, sample_bam_path) tuples
    :type in_files: list(tuple(str, str))
    :param bed_file: A path to a BED file for counting
    :type bed_file: str
    :param out_path: A path to output 
    :type out_path: [type]
    :param sample: [description], defaults to "NA"
    :type sample: str, optional
    :return: [description]
    :rtype: [type]
    """

    out_file = os.path.join(out_path, sample + COUNT_FILE_EXT)

    if os.path.exists(out_file):
        print("Loading peak counts for {s} from {f}".format(s=sample, f=out_file))
        return pd.read_csv(out_file, sep="\t")

    _nuke_these_files = {out_file: True}
    count_dfs = []

    try:
        for s_id, s_bam in in_files:
            s_out_file = os.path.join(out_path, s_id + "_coverage.bed")

            if os.path.exists(s_out_file):
                print("Loading peaks for {s} from file {f}".format(s=s_id, f=s_out_file))
                count_df = pd.read_csv(s_out_file, sep="\t", header=None)
            else:
                _nuke_these_files[s_out_file] = True

                try:

                    print("Counting peaks for {s} in file {f}".format(s=s_id, f=s_bam))
                    count_cmd = BEDTOOLS_COUNT + ["-a", bed_file, "-b", s_bam]
                    logger.debug("[BEDTOOLS COUNT]", extra={'sname': sample, 'cmd': " ".join(count_cmd)})

                    count_proc = subprocess.run(count_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                    count_df = pd.read_csv(io.StringIO(count_proc.stdout.decode("utf-8")), sep="\t", header=None)
                    count_df.to_csv(s_out_file, sep="\t", header=False, index=False)

                except (RuntimeError, pde.EmptyDataError):
                    print(" ".join(count_cmd))
                    print(io.StringIO(count_proc.stdout.decode("utf-8")).read(1000))
                    raise 

                del _nuke_these_files[s_out_file]

            count_df.columns = COUNT_COLUMNS
            count_df.rename({'count': s_id}, axis=1, inplace=True)
            count_df.drop(['num_bases_covered', 'length', 'fraction_covered'], axis=1, inplace=True)
            count_df[CHROM_COL] = count_df[CHROM_COL].astype(str)

            total_reads = _call_count_aligned(s_bam, sample=s_id) - count_df[s_id].sum()
            total_reads = pd.DataFrame([["OTHER", 1, 1, total_reads]], columns=COUNT_COLUMNS[0:3] + [s_id])

            # Put unmapped reads on the end of the 
            count_df = count_df.append(total_reads, ignore_index=True)

            count_dfs.append(count_df)

        print("Combining peak counts for {s} into file {f}".format(s=sample, f=out_file))
        full_count_df = count_dfs.pop(0)

        for cdf in count_dfs:
            full_count_df = full_count_df.merge(cdf, how='outer', on=MULTIINTER_BED_COLS).fillna(0)

        _int_cols = full_count_df.columns != CHROM_COL
        full_count_df.loc[:, _int_cols] = full_count_df.loc[:, _int_cols].astype(int)
        full_count_df.to_csv(out_file, sep="\t", index=False)
        del _nuke_these_files[out_file]

        return full_count_df

    finally:
        for k in _nuke_these_files.keys():
            _remove_file(k)
