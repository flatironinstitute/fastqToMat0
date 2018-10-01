from fastqTomat0.lib import gff
from fastqTomat0.lib import steinmetz
import argparse

ap = argparse.ArgumentParser(description="Modify GFF file transcript start and stop locations")
ap.add_argument("-g", "--gff", dest="gff_file", help="GFF FILE", metavar="FILE", required=True)
ap.add_argument("-s", "--steinmetz", dest="stein_file", help="Steinmetz 2013 FILE", metavar="FILE", required=True)
ap.add_argument("-o", "--out", dest="out_file", help="Output FILE", metavar="FILE", required=True)
args = ap.parse_args()

print("Loading new transcript data")
with open(args.stein_file, mode="r") as sf:
    stein_df = steinmetz.load_steinmetz(sf)

print("Loading GFF")
with open(args.gff_file) as gf:
    gff_df, gff_head = gff.load_gff(gf)

print("Processing new transcripts into (start,stop)")
stein_transcripts = steinmetz.process_steinmetz_data(stein_df)

print("Modifying GFF data")
gff.modify_gff_locations(gff_df, stein_transcripts)

print("Writing GFF data")
with open(args.out_file, mode="w") as of:
    gff.write_gff(of, gff_df, gff_head)


