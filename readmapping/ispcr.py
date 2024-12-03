#!/usr/bin/env python3

import subprocess

PCT_MATCH = 80.00

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    """
    main function for calling functions to perform isPCR in three steps:
    1. Identify locations where primers would anneal to the target sequence
    2. Identify pairs of locations where two primers anneal close enough together and in the correct orientation for amplification to occur
    3. Extract the amplified sequence

    Args:
        primer_file: path to the primer file
        assembly_file: path to the assembly file
        max_amplicon_size: maximum amplicon required

    Returns:
        amplicons that gets amplified in isPCR
    """
    sorted_good_hits = step_one(primer_file=primer_file, assembly_file=assembly_file)
    paired_hits = step_two(sorted_hits=sorted_good_hits, max_amplicon_size=max_amplicon_size)
    amplicons = step_three(hit_pairs=paired_hits, assembly_file=assembly_file)
    
    return amplicons

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    """
    to identify locations where primers would anneal to the target sequence

    Args:
        primer_file: file path of the primer file
        assembly_file: file path the assembly file

    Returns:
        list containing the filtered blast outputs
    """
    blast_output = call_blast(primer_file, assembly_file)
    filtered_blast_output = filter_blast(blast_output)
    
    return filtered_blast_output

def step_two(sorted_hits: list[list[str]], max_amplicon_size: int) -> list[tuple[list[str]]]:
    """
    identify pairs of locations where two primers anneal close enough together and in the correct 
    orientation for amplification to occur

    Args:
        sorted_hits: list of sorted and filtered blast outputs
        max_amplicon_size: maximum desired amplicon size

    Returns:
        a list of tuples with blast hit pairs which satisfy the condition of orientation and size
    """
    paired_hits = []
    paired_hits = find_amplicon_pairs(sorted_hits, max_amplicon_size, paired_hits)
    
    return paired_hits

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    """
    extracting amplicon sequences

    Args:
        hit_pairs: list of tuples with blast hit pairs which satisfy the condition of orientation and size
        assembly_file: file path for assembly file

    Returns:
        string with extracted amplicon sequences
    """
    bed_content = create_bed_file(hit_pairs)
    amplicon = call_seqtk(assembly_file, bed_content)
    # amplicon_formatted = pretty_print(amplicon)
    
    return amplicon[:-1] #skip the newline at the end

def call_blast(primer_file: str, assembly_file: str) -> str:
    """
    to blast primer sequence against assembly file to find sequence matches

    Args:
        primer_file: file path of the primer file
        assembly_file: file path the assembly file

    Returns:
        output of blast
    """
    blast_output = subprocess.run(["blastn", "-task", "blastn-short", "-query", primer_file, "-subject", assembly_file, "-outfmt", '6 std qlen', "-word_size=6", "-penalty=-2"], \
                   capture_output=True, \
                   text=True
                   )
    return blast_output.stdout

def filter_blast(blast_output:str) -> list[list[str]]:
    """
    filter blast hits to extracts hits which match atleast a predefined threshold PCT_MATCH and store it as a list of strings

    Args:
        blast_output: output of blast

    Returns:
        list of sorted list of blast hits with match percent above PCT_MATCH
    """
    #$3 - percent macth identity
    #$4 - match length
    #$13 - query length
    filtered_blast_output = subprocess.run("awk '{if ($3 >= PCT_MATCH && $4 == $13) print $0;}' | sort -k 9,10n", \
                                           capture_output=True, \
                                           text=True, \
                                           shell=True, \
                                           input=blast_output
                                           )
    blast_output_list = filtered_blast_output.stdout.split('\n')
    blast_output_fields = [i.split('\t') for i in blast_output_list[:-1]] #exclude last entry to get rid of unwanted newline
    return blast_output_fields

def find_amplicon_pairs(sorted_hits: list[list[str]], max_amplicon_size: int, paired_hits: list) -> list[tuple[list[str]]]:
    """
    Loop through all the sorted hits to check if any pair of hits satisfy conditions to make an amplicon
    1. Both primers anneal pointing towards one another
    2. Primers are sufficiently close to each other, set by max_amplicon_size

    Args:
        sorted_hits: list of sorted and filtered blast outputs
        max_amplicon_size (int): maximum desired amplicon size
        paired_hits (list): list to store hit pairs

    Returns:
        list of tuples with blast hit pairs which satisfy the condition of orientation and size
    """
    #primer[1] - query seqID
    #primer[0] - subject seqID
    #primer[8] and [9] - matched start and end position
    for primer1 in sorted_hits:
        for primer2 in sorted_hits:
            valid_amplicon_pair = ()            
            # We can skip comparing the same hits as they cannot make an amplicon
            if primer1 == primer2:
                continue
            else:
                # Check if both hits have the same sequence ID and primer IDs are not the same
                if primer1[1] == primer2[1] and primer1[0] != primer2[0]: 
                    # Compare the 3' end of both primers and if their difference is less than amplicon size, they are valid pairs
                    if int(primer1[9]) < int(primer2[9]) and int(primer2[9]) - int(primer1[9]) < max_amplicon_size:
                        valid_amplicon_pair = (primer1, primer2)
                        paired_hits.append(valid_amplicon_pair)

    return paired_hits

def create_bed_file(hit_pairs: list[tuple[list[str]]]) -> str:
    """
    create bed file using the filtered amplicon list

    Args:
        hit_pairs: list containing hit pairs

    Returns:
        BED content as string
    """
    bed_list = []
    #extract seqID, amplicon co-ordinate which is the 3' end position of the two primers
    for amplicon_pair in hit_pairs:
        primer1, primer2 = amplicon_pair
        bed_list.append(f"{primer1[1]}\t{primer1[9]}\t{int(primer2[9])-1}")
    
    #create a string containing BED contents
    bed_content = ('\n').join(bed_list)
    return bed_content

def call_seqtk(assembly_file: str, bed_content: str) -> str:
    """
    execute seqtk to extract sequence from assembly file using coordinates in bed file

    Args:
        assembly_file: path to assembly file
        bed_file: BED content

    Returns:
        string of amplicon sequences extracted from hit positions
    """
    amplicon_seq = subprocess.run(f'seqtk subseq {assembly_file} <(echo "{bed_content}" ; data/Vibrio_cholerae_N16961.bed | xargs)', \
                   shell=True, \
                   capture_output=True, \
                   text=True, \
                   executable="/bin/bash"
                   )
    return amplicon_seq.stdout

def pretty_print(string: str, position: int = 83) -> str:
    """
    Break the sequence into readable lengths by adding newline character

    Args:
        string : long string to be broken
        position : index at which new line character should be added. Defaults to 83.
    
    Returns:
        String with added new line characters
    """
    str_list = string.split('\n')
    #for every even element in the list containing sequence, break the sequence into 'position' number of bases
    for line in range(1, len(str_list), 2):
        formatted_str = '\n'.join(str_list[line][i:i+position] for i in range(0, len(str_list[line]), position))
        str_list[line] = formatted_str
    
    return '\n'.join(str_list)