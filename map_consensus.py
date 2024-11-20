#!/usr/bin/env python3

import argparse
import subprocess
from readmapping import sam_record_processing

#accept and parse command-line arguments
parser = argparse.ArgumentParser(prog="map_consensus.py", 
                                 description="Map reads to reference to get the best mapping and extract consensus sequence")

read1_file = parser.add_argument(
    "-1", "--read1",
    help="First paired-end FASTQ file path",
    dest="read1",
    type=str,
    required=True
)

read2_file = parser.add_argument(
    "-2", "--read2",
    help="Second paired-end FASTQ file path",
    dest="read2",
    type=str,
    required=True
)

ref_file = parser.add_argument(
    "-r", "--reference",
    help="Reference FASTA file path",
    dest="ref",
    type=str,
    required=True
)

seq_name_arg = parser.add_argument(
    "-s", "--reference-name",
    help="[Optional] Sequence name in reference",
    dest="seq_name",
    type=str,
    required=False,
)

args = parser.parse_args()

#map reads to reference
cmd_list = ["minimap2", "-ax", "sr", "-B", "0", "-k", "10", args.ref, args.read1, args.read2]
map_output = subprocess.run(cmd_list, capture_output=True, text=True)

if map_output.returncode == 0:
    accession = args.read1.split("/")[-1].split("_")[0]
    with open(f"ex11_data/{accession}.sam", "w") as sam_writer:
        sam_writer.write(map_output.stdout)
else:
    raise Exception("Mapping of reads returned non-zero exit code")

sam_path = f"ex11_data/{accession}.sam"
sam = sam_record_processing.SAM.from_sam(sam_path)
if args.seq_name: #if optional sequence name is provided
    consensus_str = sam.consensus(args.seq_name)
    print(f">{args.seq_name}_consensus\n{consensus_str}")
else: #print best consensus
    consensus_str = sam.best_consensus()
    print(f"{consensus_str}")