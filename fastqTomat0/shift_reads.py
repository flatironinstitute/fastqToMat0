import argparse
import glob
import os
import pathlib
import multiprocessing
import itertools
import gzip
import subprocess


def main():
    ap = argparse.ArgumentParser(description="Shift UMI sequences together for downstream processing")
    ap.add_argument("-f", "--fastq", dest="fastq", help="Paired-End Read Barcode FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-o", "--out", dest="out", help="Output path", metavar="PATH", required=True)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    ap.add_argument("--pattern", dest="pattern", help="Barcode & UMI Pattern (BBBBUUU)", metavar="PATTERN",
                    default="UUUUBBBBBBBBBBBBUUUUUUUUUU")
    args = ap.parse_args()

    file_names = []
    for fastq in args.fastq:
        file_names.extend(glob.glob(fastq))

    for fn in file_names:
        if not os.path.exists(fn):
            raise FileNotFoundError("{fq} File Not Found".format(fq=fn))

    out_path = pathlib.Path(os.path.expanduser(args.out))
    out_path.mkdir(parents=True, exist_ok=True)

    pattern = args.pattern.upper()

    if args.cores == 1:
        list(map(shift_reads, file_names, itertools.repeat(out_path), itertools.repeat(pattern)))
    else:
        with multiprocessing.Pool(processes=args.cores) as mp:
            mp.starmap(shift_reads, zip(file_names, itertools.repeat(out_path), itertools.repeat(pattern)))


def shift_reads(file, out_path, pattern, gz=False):

    file = os.path.abspath(file)
    _, file_name = os.path.split(file)

    _breakpoints = [0] + [i for i in range(1, len(pattern)) if pattern[i] != pattern[i-1]] + [len(pattern)]
    _plen = len(pattern)

    bc, um = [], []
    for i in range(1, len(_breakpoints)):
        slicer = "x[{start}:{stop}]".format(start=_breakpoints[i - 1], stop=_breakpoints[i])
        bc.append(slicer) if all(c == "B" for c in pattern[_breakpoints[i - 1]: _breakpoints[i]]) else um.append(slicer)

    def _shifter(x): return eval(" + ".join(bc) + " + " + " + ".join(um))

    _chrp = list(map(chr, range(65, 65 + len(pattern))))
    print("Shifting pattern {p} to {np} ({chrp} to {chrpn})".format(p=pattern, np=_shifter(pattern),
                                                                    chrp=_chrp, chrpn=_shifter(_chrp)))

    o_file = os.path.join(out_path, file_name[:-3] if file_name.endswith(".gz") else file_name)
    infh = gzip.open(file, mode="rt", encoding="utf-8") if file_name.endswith(".gz") else open(file)

    with open(o_file, mode="wt") as outfh:
        for i, li in enumerate(infh):

            if i % int(1e6) == 0:
                print("{f}: {i} records parsed".format(f=file_name, i=i))

            li = li.strip()
            _no_process = li.startswith("@") or (len(li) == 0) or li.startswith("+")
            if _no_process:
                print(li, file=outfh)

            else:
                if len(li) != _plen:
                    _err = "Line {i}: {v} is {ll} characters, pattern is {pl}".format(i=i, v=li, ll=len(li), pl=_plen)
                    raise ValueError(_err)
                print(_shifter(li), file=outfh)

        print("{f} shift complete".format(f=file_name))

        if gz:
            print("Zipping output file {fi}".format(fi=o_file))
            subprocess.call(["gzip", os.path.abspath(os.path.expanduser(os.path.expandvars(o_file)))])

    infh.close()


if __name__ == '__main__':
    main()
