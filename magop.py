#!/usr/bin/env python3

import argparse
import subprocess
import os
import tempfile
import re
from readmapping import sam_record_processing, ispcr, nw

def reverse_complement(string:str) -> str:
    """
    Generate the reverse complement of the string using complement_dict

    Args:
        string : string to be reverse complemented

    Returns:
        reverse complemented string
    """
    complement_dict = {
    "A":"T",
    "G":"C",
    "T":"A",
    "C":"G",
    "N":"N"
    }
    return "".join(complement_dict.get(base) for base in string[::-1])

#accept and parse command-line arguments
parser = argparse.ArgumentParser(prog="magop.py", 
                                 description=f'''Infer phylogenetic relationship among given organisms using paired-end reads or assemblies\n
                                        Usage: magop.py [-a ASSEMBLY [ASSEMBLY ...]] -p PRIMERS [-r READS [READS ...]] [-s REF_SEQS]''')

assembly_files = parser.add_argument(
    "-a", "--assemblies",
    help="Assembly files as space separated list or directory conatining them in FASTA format",
    dest="assemblies",
    type=str,
    nargs='*',
    required=False
)

read_files = parser.add_argument(
    "-r", "--reads",
    help="PE Illumina reads as space separated list or directory conatining them in FASTQ format",
    dest="reads",
    type=str,
    nargs='*',
    required=False
)

primer_file = parser.add_argument(
    "-p", "--primer-file",
    help="Path to the primer file",
    dest="primer",
    type=str,
    nargs=1,
    required=True
)

ref_file = parser.add_argument(
    "-s", "--reference",
    help="File path to Reference FASTA. Mandatory with FASTQ file inputs",
    dest="ref",
    type=str,
    nargs='*',
    required=False
)


args = parser.parse_args()
amplicon_list = []

# print(args.assemblies, args.reads, args.primer, args.ref)
if args.assemblies:
    for assembly_file in args.assemblies:
        # print(assembly_file)
        amplicons = ispcr.ispcr(args.primer[0], assembly_file, 2000) 
        asm_amp_seq = amplicons.split('\n')[-1] #get one amplicon
        asm_amp_header = amplicons.split('\n')[-2].split(':')[0]
        amplicon_list.append(f"{asm_amp_header}\n{asm_amp_seq}")
    # print(amplicon_list)
        # break


if args.reads:
    if not args.ref:
        parser.error("Reference file is mandatory with FASTQ input")
    matching_reads = {}
    for read_file in args.reads:
        accession = (re.findall(r"[A-Za-z0-9]*/*([A-Za-z0-9]+)_[12].fastq", read_file))[0]

        if accession in matching_reads:
            matching_reads[accession].append(read_file)
        else:
            matching_reads[accession] = [read_file]
    # print(matching_reads)
    for accession, pe_reads in matching_reads.items():
        if len(pe_reads) != 2 or ("_1" not in pe_reads[0] and "_2" not in pe_reads[1]): #error if reads are not paired
            parser.error(f"Read files for {accession} is not paired. Input paired-end FASTQ files")
        

    ref_16S = "ex11_data/refs/16S.fna"

    for accession, pe_reads in matching_reads.items():
        # print(args.ref[0], pe_reads[0], pe_reads[1])
        cmd_list = ["minimap2", "-ax", "sr", "-B", "0", "-k", "10", args.ref[0], pe_reads[0], pe_reads[1]]
        map_output = subprocess.run(cmd_list, capture_output=True, text=True)

        if map_output.returncode == 0:
            with open(f"ex12_data/sams/{accession}.sam", "w") as sam_writer: #delete the folder and update to tempdir
                sam_writer.write(map_output.stdout)
        else:
            raise Exception("Mapping of reads returned non-zero exit code")

        sam_path = f"ex12_data/sams/{accession}.sam"
        sam = sam_record_processing.SAM.from_sam(sam_path)
        consensus_str = sam.best_consensus()
        # print(consensus_str)

        
        with tempfile.NamedTemporaryFile(mode='t+w') as f_cons:
            f_cons.write(consensus_str)
            f_cons.seek(0)
            
            #perform isPCR to extract v4 region from 16S
            amplicons = ispcr.ispcr(args.primer[0], f_cons.name, 2000) 
            amplicon_v4 = amplicons.split('\n')[-1] #get one amplicon
            amplicon_list.append(f">{accession}\n{amplicon_v4}")
            # print(f"{accession} amplicon: {amplicon_v4}")
            

    # isPCR reference v4 sequences
    ref_amplicons = ispcr.ispcr(args.primer[0], args.ref[0], 2000)

    for ref_amp in ref_amplicons.split('\n'):
        if ref_amp.startswith('>'):
            ref_header = ""
            ref_header = ref_amp.split(':')[0]
        else:
            amplicon_list.append(f"{ref_header}\n{ref_amp}")
    # print(ref_amplicons)

    # print(amplicon_list)

    aligned_amplicons = []
    for amplicon in amplicon_list:
        _, anchor_seq = amplicon_list[0].split('\n')
        amp_header, amp_seq = amplicon.split('\n')
        # print(anchor_seq, amp_header, amp_seq)
        _, score = nw.needleman_wunsch(amp_seq, anchor_seq, 1, -1, -1)
        _, rev_score = nw.needleman_wunsch(reverse_complement(amp_seq), anchor_seq, 1, -1, -1)

        if score > rev_score:
            aligned_amplicons.append(amplicon)
        else:
            aligned_amplicons.append(f"{amp_header}\n{reverse_complement(amp_seq)}")
            # print(f">{amp_header}\n{reverse_complement(amp_seq)}")

for aln_amp in aligned_amplicons:
    print(aln_amp)


# ./magop.py -a ex12_data/assemblies/* -p ex12_data/primers/general_16S_515f_806r.fna -r ex12_data/reads/* -s ex12_data/refs/V4.fna > readmapping/aligned_v4.fna   
# muscle -align readmapping/aligned_v4.fna -output readmapping/multi_aligned_v4.txt
# iqtree2 -s readmapping/multi_aligned_v4.txt --prefix readmapping/ml_tree --quiet -redo
# gotree reroot midpoint -i readmapping/ml_tree.treefile -o readmapping/rooted_ml_tree





        



