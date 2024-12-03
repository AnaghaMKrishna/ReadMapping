from __future__ import annotations

class SAM:
    def __init__(self, primary_reads:list[SAMRead], references:dict[str, list[int, int, int]]):
        """
        initialize the sam records in the SAM file

        Args:
            primary_reads: list containing SAMRead object for primary, non header lines in SAM file
            references: a dictionary of reference name, it's start and end coordinates and number of mapped reads
        """
        self.primary_reads = primary_reads
        self.references = references

    @classmethod
    def from_sam(cls, file_path:str) -> SAM:
        """
        creates an instance of SAM class from SAM file

        Args:
            file_path: SAM file path

        Returns:
            an instance of the SAM class 
        """
        primary_reads = []
        ref_dict = {}
        with open(file_path) as sam_file:
            for line in sam_file:
                if not line.startswith('@'): 
                    sam_line = SAMRead(line.strip())
                    if sam_line.is_primary and sam_line.is_mapped:
                        primary_reads.append(sam_line) #store only primary reads
                        ref_dict[sam_line.rname][2] += 1 #increment count for read mapped to a reference
                elif line.startswith('@SQ'): #@SQ	SN:ref_seq	LN:ref_len
                    header_line = line.strip().split('\t')
                    ref_seq = header_line[1] #value of SN:
                    ref_len = header_line[2] #value of LN:
                    ref_dict[ref_seq[3:]] = [1, int(ref_len[3:]), 0] #skip SN: and LN:
                else:
                    continue
        sam = cls(primary_reads, ref_dict) 
        return sam

    def reads_at_pos(self, rname:str, pos:int)-> list[SAMRead]:
        """
        finds all reads mapped at a position on reference

        Args:
            rname: reference sequence name
            pos: position relative to the reference

        Returns:
            a list of read objects mapped at the position
        """
        reads_list = []
        for read in self.primary_reads:
            if read.rname == rname and read.start <= pos <= read.end:
                reads_list.append(read)

        return reads_list
    
    def pileup_at_pos(self, rname:str, pos:int)-> tuple[list[str], list[str]]:
        """
        finds all bases and corresponding quality score mapped at a position on reference

        Args:
            rname: reference sequence name
            pos: position relative to the reference

        Returns:
            a tuple of list of mapped bases and list of corresponding quality scores
        """
        base_list = []
        qual_list = []
        for read in self.reads_at_pos(rname, pos):
            base_list.append(read.base_at_pos(pos))
            qual_list.append(read.qual_at_pos(pos))

        return base_list, qual_list

    def consensus_at_pos(self, rname:str, pos:int) -> str:
        """
        finds the consensus base at a position on the reference

        Args:
            rname: reference sequence name
            pos: position relative to the reference

        Returns:
            a base with highest count at the position
        """
        
        bases_at_pos = self.pileup_at_pos(rname, pos)[0]
        base_count_dict = {base: bases_at_pos.count(base) for base in bases_at_pos}
        base_counts = sorted(base_count_dict.items(), key=lambda item: item[1]) #sort based on counts

        if len(base_counts) > 1: #multiple bases
            for base, count in base_counts:
                if not count == base_counts[-1][1]: #check if count is equal to max_count
                    continue
                if base_counts[-1][1] == base_counts[-2][1]: #tie between bases
                    return "N"
                if count > 0.5 * len(bases_at_pos): #call base with >50% majority
                    return base
                else: #no majority base
                    return "N"
        elif len(base_counts) == 1: #all reads agree with the base
            return base_counts[0][0]

        return "" #no reads mapped or rname does not exist
    
    def consensus(self, rname:str) -> str:
        """
        consensus sequence for a reference sequence

        Args:
            rname: reference sequence name

        Returns:
            consensus sequence of all positions in the reference
        """
        if rname not in self.references: #if reference name does not exist
            return ""
        cons_base_list = []
        ref_start = self.references[rname][0]
        ref_end = self.references[rname][1]
        for i in range(ref_start, ref_end + 1): #1-base indexing
            cons_base = self.consensus_at_pos(rname, i)
            cons_base_list.append(cons_base)
        return "".join(cons_base_list) if cons_base_list else "" # empty string if no reads are mapped to rname
    
    def best_consensus(self) -> str:
        """
        finds the best consensus across all reference sequences

        Returns:
            best consensus sequence
        """
        #find reference with maximum number of mapped reads
        best_rname = max(self.references, key=lambda ref: self.references[ref][2])
        best_cons = self.consensus(best_rname)

        return f">{best_rname}_consensus\n{best_cons}"


class SAMRead:
    def __init__(self, record:str):
        """
        initialize the fields in record and define attributes

        Args:
            record: aligned read in SAM file
        """
        record = record.strip().split("\t")
        self.qname: str = record[0] 
        self.flag:int = int(record[1])
        self.rname:str = record[2] 
        self.pos:int = int(record[3])
        self.mapq:str = record[4] 
        self.cigar:str = record[5] 
        self.rnext:str = record[6] 
        self.pnext:int = int(record[7])
        self.tlen:int = int(record[8])
        self.seq:str = record[9] 
        self.qual:str = record[10]
        self.extras:list[str] = record[11:]

        #decode flag attributes
        self.bin_flag = f"{self.flag:012b}"[::-1] #reverse to get LSB at index 0
        self.is_mapped = True if self.bin_flag[2] == '0' else False
        self.is_forward = True if self.bin_flag[4] == '0' else False
        self.is_reverse = False if self.bin_flag[4] == '0' else True
        self.is_primary = True if self.bin_flag[8] == '0' and self.bin_flag[11] == '0' else False
        
        #list of mapped bases and quality score of bases in mapped reads
        self.mapped_bases, self.mapped_qual = self.find_mapped_elements() if self.is_mapped else ([], [])

        #calculate read range
        self.start = self.pos
        self.end = self.start + len(self.mapped_bases) - 1
        
    def process_cigar(self) -> tuple[list[list[str, int]], int]:
        """
        read cigar into a list and offset clipped bases

        Returns:
            cigar list and offset index if read is soft clipped
        """
        num = ''
        prev_num = 0
        cigar_list = []
        clip_idx = 0

        #convert cigar string into list - [[character, number]]
        for c in self.cigar:
            if c not in 'MDISH':
                num += c
            else: 
                cigar_list.append([c, prev_num + int(num)])
                prev_num = cigar_list[-1][1]
                num = ''
        # print(cigar_list)
        
        #update cigar list after offsetting soft clipped bases from start or end of read
        if cigar_list[0][0] == 'S':
            clip_idx = cigar_list[0][1]
            cigar_list = [[l[0], l[1] - clip_idx] for l in cigar_list]
            return cigar_list, clip_idx
        if cigar_list[-1][0] == 'S':
            clip_idx = cigar_list[-2][1] 
            cigar_list = [[l[0], l[1] - clip_idx] for l in cigar_list]
            return cigar_list, clip_idx 
        
        return cigar_list, 0

    def find_mapped_elements(self) -> tuple[list[str], list[str]]:
        """
        returns a list of mapped elements relative to the reference

        Returns:
            a tuple of list containing mapped bases and mapped quality scores
        """

        base_list = [base for base in self.seq]
        qual_list = [qual for qual in self.qual]
        cigar_list, clip_idx = self.process_cigar()

        #remove soft clipped elements
        if cigar_list[0][0] == 'S': 
            del base_list[:clip_idx]
            del qual_list[:clip_idx]
        if cigar_list[-1][0] == 'S':
            del base_list[clip_idx:]
            del qual_list[clip_idx:]

        to_del = []
        #process the element list according to cigar
        for n,l in enumerate(cigar_list):
            if l[0] == 'S':
                continue #already handled

            elif l[0] == 'H':
                continue #no action required
            
            elif l[0] == 'M':
                continue #no action required
            
            elif l[0] == 'D': #add empty string to denote deleted bases - [ele1, ele2, ' ', ele3]
                start_idx = cigar_list[n-1][1]
                stop_idx = l[1]
                diff = stop_idx - start_idx
                base_list[start_idx:start_idx] = list(' ' * diff)
                qual_list[start_idx:start_idx] = list(' ' * diff)

            elif l[0] == 'I': #combine bases into list to denote insertion - [ele1, [ele2, ele3], ele4]
                start_idx = cigar_list[n-1][1] - 1
                stop_idx = l[1] 
                diff = stop_idx - start_idx - 1
                base_list[start_idx] = base_list[start_idx:stop_idx]
                qual_list[start_idx] = qual_list[start_idx:stop_idx]
                to_del.append([start_idx + 1, diff]) #keep track of index and length to remove bases after combining into list

        #remove elements next to insertions
        for ele in to_del[::-1]: #iterate from end to avoid changing list index
            del base_list[ele[0]:ele[0]+ele[1]]
            del qual_list[ele[0]:ele[0]+ele[1]]

        #get bases and quality scores for indices present in idx_list
        self.mapped_bases = base_list
        self.mapped_qual = qual_list

        return self.mapped_bases, self.mapped_qual

    def base_at_pos(self, pos: int) -> str:
        """
        return the base mapped at a position

        Args:
            pos: position relative to the reference

        Returns:
            base at pos
        """
        base = ""

        # unmapped read
        if not self.is_mapped:
            return ""

        #check if pos within the range of mapped read
        if self.start <= pos <= self.end:
            idx = abs(self.pos - pos) 
            base = "".join(self.mapped_bases[idx])
            return "" if base == ' ' else base #bases at deletion
        else:
            return ""

    def qual_at_pos(self, pos: int) -> str:
        """
        return the quality score mapped at a position

        Args:
            pos: position relative to the reference

        Returns:
            quality at pos
        """
        qual = ""

        #unmapped read
        if not self.is_mapped:
            return ""

        #check if pos within the range of mapped read
        if self.start <= pos <= self.end:
            idx = abs(self.pos - pos) 
            qual = "".join(self.mapped_qual[idx])
            return "" if qual == ' ' else qual
        else:
            return ""

    def mapped_seq(self) -> str:
        """
        returns mapped sequence of the read

        Returns:
            sequence mapped to reference
        """
        base_list = self.mapped_bases

        if 'I' in self.cigar:
            base_list = ["".join(base) for base in base_list] #flatten 2D list with insertions into 1D
        elif 'D' in self.cigar:
            base_list = ["-" if base == " " else base for base in base_list] #replace deletions " " with "-"
        return "".join(base_list)
