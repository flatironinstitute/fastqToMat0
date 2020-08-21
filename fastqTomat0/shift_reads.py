import argparse
import glob
import os
import pathlib
import multiprocessing
import itertools
import gzip
import subprocess


def main():
    ap = argparse.ArgumentParser(description="Coerce non-10x barcodes to look like 10x barcodes (to test cellranger)")
    ap.add_argument("-f", "--fastq", dest="fastq", help="Paired-End Read Barcode FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-o", "--out", dest="out", help="Output path", metavar="PATH", required=True)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    args = ap.parse_args()

    file_names = []
    for fastq in args.fastq:
        file_names.extend(glob.glob(fastq))

    for fn in file_names:
        if not os.path.exists(fn):
            raise FileNotFoundError("{fq} File Not Found".format(fq=fn))

    out_path = pathlib.Path(os.path.expanduser(args.out))
    out_path.mkdir(parents=True, exist_ok=True)

    if args.cores == 1:
        list(map(shift_reads, file_names, itertools.repeat(out_path)))
    else:
        with multiprocessing.Pool(processes=args.cores) as mp:
            mp.starmap(shift_reads, zip(file_names, itertools.repeat(out_path)))


def shift_reads(file, out_path, shift_n=4, shift_loc=12, gz=True):

    file = os.path.abspath(file)
    _, file_name = os.path.split(file)

    a, b = shift_n, shift_n + shift_loc

    l_idx = 0

    o_file = os.path.join(out_path, file_name[:-3] if file_name.endswith(".gz") else file_name)
    infh = gzip.open(file, mode="rt", encoding="utf-8") if file_name.endswith(".gz") else open(file)

    with open(o_file, mode="wt") as outfh:
        for i, li in enumerate(infh):

            if i % int(1e6) == 0:
                print("{f}: {i} records parsed".format(f=file_name, i=i))

            li = li.strip()
            _no_process = li.startswith("@") or (len(li) == 0) or li.startswith("+")
            print(li, file=outfh) if _no_process else print(li[a:b] + li[0:a] + li[b:], file=outfh)

        print("{f} shift complete".format(f=file_name))

        if gz:
            print("Zipping output file {fi}".format(fi=o_file))
            subprocess.call(["gzip", os.path.abspath(os.path.expanduser(os.path.expandvars(o_file)))])

    infh.close()


if __name__ == '__main__':
    main()
