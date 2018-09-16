from fastqTomat0.fastqBCLinker import link_barcodes

is_zipped = True

BC1 = ["/scratch/cj59/10x_out/Undetermined_S0_L001_R1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L002_R1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L003_R1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L004_R1_001.fastq.gz"] # 26bp cell index
BC2 = ["/scratch/cj59/10x_out/Undetermined_S0_L001_R2_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L002_R2_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L003_R2_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L004_R2_001.fastq.gz"] # 98bp transcript index
BC3 = ["/scratch/cj59/10x_out/Undetermined_S0_L001_I1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L002_I1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L003_I1_001.fastq.gz",
       "/scratch/cj59/10x_out/Undetermined_S0_L004_I1_001.fastq.gz"] # 8bp library index

transcript_barcdode_fasta = "/home/cj59/Documents/all_barcodes_trim.fa"

indexes = ["ATTACTCG",
           "TCCGGAGA"] # Index sequences to retain

max_index_mismatches = 1

output_file_path = "/scratch/cj59/10x_out/bc_decode.tsv"

bc1_pattern = "BBBBBBBBBBBBBBBBUUUUUUUUUU"
bc2_seed = "ttttctaaggatcca"
bc2_len = 19

link_barcodes(BC1, BC2, BC3, output_file_path, is_zipped=is_zipped, bc1_pattern=bc1_pattern, bc2_seed=bc2_seed,
              bc2_len=bc2_len, allowed_indexes=indexes, max_index_mismatch=max_index_mismatches,
              bc2_mapfile=transcript_barcdode_fasta)