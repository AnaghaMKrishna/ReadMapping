#!/usr/bin/env python3

import argparse
import subprocess
from collections import defaultdict
from tempfile import TemporaryDirectory, NamedTemporaryFile
import os
import re
from readmapping import process_sam, ispcr, nw

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

def handle_assemblies(primer_path:str, assemblies_path:list[str], max_amplicon_size:int) -> list[str]:
    """
    perform isPCR to extract amplicons from assembly files

    Args:
        primer_path: path to FASTA formatted primer file
        assemblies_path: path to FASTA formatted assembly files

    Returns:
        list of amplicons extracted from assemblies
    """
    assembly_amplicons = []
    for assembly_file in assemblies_path:
        amplicons = ispcr.ispcr(primer_path, assembly_file, max_amplicon_size) 
        asm_amp_seq = amplicons.split('\n')[-1] #get one amplicon
        asm_amp_header = amplicons.split('\n')[-2].split(':')[0]
        assembly_amplicons.append(f"{asm_amp_header}\n{asm_amp_seq}")
    
    return assembly_amplicons

def handle_reads(primer_path:str, reads_path:list[str], ref_path:str, max_amplicon_size:int) -> list[str]:
    """
    perform isPCR on consensus sequence obtained from paired-end read files

    Args:
        primer_path: path to FASTA formatted primer file
        reads_path: path to FASTQ formatted read files

    Returns:
        list of amplicons extracted from read files
    """
    read_amplicons = []
    matching_reads = defaultdict(list)

    #extract accession and match paired read files
    for read_file in reads_path:
        accession = (re.findall(r"[A-Za-z0-9]*/*([A-Za-z0-9]+)_[12].fastq", read_file))[0]
        matching_reads[accession].append(read_file)
    for accession, pe_reads in matching_reads.items():
        if len(pe_reads) != 2 or ("_1" not in pe_reads[0] and "_2" not in pe_reads[1]): #error if reads are not paired
            raise Exception(f"Read files for {accession} is not paired. Input paired-end FASTQ files")        

    with TemporaryDirectory() as sam_tempdir:
        for accession, pe_reads in matching_reads.items():
            cmd_list = ["minimap2", "-ax", "sr", "-B", "0", "-k", "10", ref_path, pe_reads[0], pe_reads[1]]
            map_output = subprocess.run(cmd_list, capture_output=True, text=True)

            #save SAM files in temporary directory
            sam_path = os.path.join(sam_tempdir, f"{accession}.sam")
            if map_output.returncode == 0:
                with open(sam_path, "w") as sam_writer:
                    sam_writer.write(map_output.stdout)
                    sam = process_sam.SAM.from_sam(sam_path)
                    consensus_str = sam.best_consensus()

                    #write consensus sequence to temporary file
                    with NamedTemporaryFile(mode='t+w') as f_cons:
                        f_cons.write(consensus_str)
                        f_cons.seek(0)
                        amplicons = ispcr.ispcr(primer_path, f_cons.name, max_amplicon_size) 
                        amplicon_v4 = amplicons.split('\n')[-1] #get one amplicon
                        read_amplicons.append(f">{accession}\n{amplicon_v4}")

            else:
                raise Exception(f"Mapping of reads returned non-zero exit code:\n{map_output.stderr}")    
    
    return read_amplicons

def handle_references(primer_path:str, ref_path:list[str], max_amplicon_size:int) -> list[str]:
    """
    perform isPCR to extract amplicons from reference file

    Args:
        primer_path: path to FASTA formatted primer file
        ref_path: path to multi-FASTA reference file

    Returns:
        list of amplicons extracted from reference sequences
    """
    ref_amplicons = []
    amplicons = ispcr.ispcr(primer_path, ref_path, max_amplicon_size)

    for ref_amp in amplicons.split('\n'):
        if ref_amp.startswith('>'):
            ref_header = ""
            ref_header = ref_amp.split(':')[0]
        else:
            ref_amplicons.append(f"{ref_header}\n{ref_amp}")
    
    return ref_amplicons

def orient_amplicons(amplicon_list: list[str], match: int, mismatch: int, gap: int)-> list[str]:
    """
    orient all amplicons in same direction

    Args:
        amplicon_list: list containing extracted amplicons
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        list of oriented amplicons
    """

    oriented_amplicons = []
    for amplicon in amplicon_list:
        _, anchor_seq = amplicon_list[0].split('\n')
        amp_header, amp_seq = amplicon.split('\n')

        _, score = nw.needleman_wunsch(amp_seq, anchor_seq, match, mismatch, gap)
        _, rev_score = nw.needleman_wunsch(reverse_complement(amp_seq), anchor_seq, match, mismatch, gap)

        if score > rev_score:
            oriented_amplicons.append(amplicon)
        else:
            oriented_amplicons.append(f"{amp_header}\n{reverse_complement(amp_seq)}")

    return oriented_amplicons

def process_args():
    """
    parse command-line arguments

    Returns:
        Namespace containing parsed arguments
    """
    #accept and parse command-line arguments
    parser = argparse.ArgumentParser(prog="magop.py", 
                                    description=f'''Infer phylogenetic relationship among given organisms using paired-end reads or assemblies\n
                                            Usage: magop.py [-a ASSEMBLY [ASSEMBLY ...]] -p PRIMERS [-r READS [READS ...]] [-s REF_SEQS]
                                            [--max-amplicon-size] [--match] [--mismatch] [--gap]''')

    parser.add_argument(
        "-a", "--assemblies",
        help="Assembly files as space separated list or directory conatining them in FASTA format",
        dest="assemblies",
        type=str,
        nargs='*',
        required=False
    )

    parser.add_argument(
        "-r", "--reads",
        help="PE Illumina reads as space separated list or directory conatining them in FASTQ format",
        dest="reads",
        type=str,
        nargs='*',
        required=False
    )

    parser.add_argument(
        "-p", "--primer-file",
        help="Path to the primer file",
        dest="primer",
        type=str,
        nargs=1,
        required=True
    )

    parser.add_argument(
        "-s", "--reference",
        help="File path to Reference FASTA. Mandatory with FASTQ file inputs",
        dest="ref",
        type=str,
        nargs='*',
        required=False
    )

    parser.add_argument(
        "--max-amplicon-size",
        help="maximum amplicon size for isPCR",
        dest="max_amplicon_size",
        type=int,
        default=2000,
        required=False
    )

    parser.add_argument(
        "--match",
        help="match score to use in alignment",
        dest="match",
        type=int,
        default=1,
        required=False
        )

    parser.add_argument(
        "--mismatch",
        help="mismatch penalty to use in alignment",
        dest="mismatch",
        type=int,
        default=-1,
        required=False
    )

    parser.add_argument(
        "--gap",
        help="gap penalty to use in alignment",
        dest="gap",
        type=int,
        default=-1,
        required=False
    )

    args = parser.parse_args()
    return args

def main():
    args = process_args()
    amplicon_list = []

    if args.assemblies:
        assembly_amplicons = handle_assemblies(args.primer[0], args.assemblies, args.max_amplicon_size)
        amplicon_list.extend(assembly_amplicons)

    if args.reads:
        if not args.ref:
            raise Exception("Reference file is mandatory with FASTQ input")
        reads_amplicons = handle_reads(args.primer[0], args.reads, args.ref[0], args.max_amplicon_size)
        amplicon_list.extend(reads_amplicons)
        ref_amplicons = handle_references(args.primer[0], args.ref[0], args.max_amplicon_size)
        amplicon_list.extend(ref_amplicons)

    oriented_amplicons = orient_amplicons(amplicon_list, args.match, args.mismatch, args.gap)
    for amp in oriented_amplicons:
        print(amp)

if __name__ == "__main__":
    main()

#USAGE:
# ./phyl_analysis.py -a ex12_data/assemblies/* -p ex12_data/primers/general_16S_515f_806r.fna -r ex12_data/reads/* -s ex12_data/refs/V4.fna > readmapping/outputs/aligned_v4.fna   
# multi-sequence alignment  - muscle -align readmapping/outputs/aligned_v4.fna -output readmapping/outputs/multi_aligned_v4.txt
# generate phylogenetic tree - iqtree2 -s readmapping/outputs/multi_aligned_v4.txt --prefix readmapping/outputs/ml_tree --quiet -redo
# mid-point rooting tree - gotree reroot midpoint -i readmapping/outputs/ml_tree.treefile -o readmapping/outputs/rooted_ml_tree