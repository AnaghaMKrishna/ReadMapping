#!/usr/bin/env python3

from __future__ import annotations

class SAM:
    def __init__(self, primary_reads:list[SAMRead]):
        """
        initialize the sam records in the SAM file

        Args:
            file_path: SAM file path
        """
        self.primary_reads = primary_reads

    @classmethod
    def from_sam(cls, file_path:str) -> SAM:
        primary_reads = []
        with open(file_path) as sam_file:
            for line in sam_file:
                if not line.startswith('@'): 
                    sam_line = SAMRead(line.strip())
                    primary_reads.append(sam_line) if sam_line.is_primary else '' #store only primary reads
        # print(primary_reads)
        sam = cls(primary_reads)  
        return sam

    def reads_at_pos(self, rname:str, pos:int)-> list[SAMRead]:
        reads_list = []
        for read in self.primary_reads:
            if read.rname == rname and read.read_range[0] <= pos <= read.read_range[1]:
                reads_list.append(read)

        # print(reads_pos_list)
        return reads_list
    
    def pileup_at_pos(self, rname:str, pos:int)-> tuple[list[str], list[str]]:
        base_list = []
        qual_list = []
        for read in self.reads_at_pos(rname, pos):
            base_list.append(read.base_at_pos(pos))
            qual_list.append(read.qual_at_pos(pos))

        # print(base_list, qual_list)
        return base_list, qual_list


class SAMRead:
    def __init__(self, record:str):
        """
        initialize the fields in record and define attributes

        Args:
            record: aligned read in SAM file
        """
        record = record.strip().split("\t")
        self.qname = record[0] 
        self.flag = int(record[1])
        self.rname = record[2] 
        self.pos = int(record[3])
        self.mapq = record[4] 
        self.cigar = record[5] 
        self.rnext = record[6] 
        self.pnext = int(record[7])
        self.tlen = int(record[8])
        self.seq = record[9] 
        self.qual = record[10]
        self.extras = record[11:]

        #decode flag attributes
        self.bin_flag = f"{self.flag:012b}"[::-1] #reverse to get LSB at index 0
        self.is_mapped = True if self.bin_flag[2] == '0' else False
        self.is_forward = True if self.bin_flag[4] == '0' else False
        self.is_reverse = False if self.bin_flag[4] == '0' else True
        self.is_primary = True if self.bin_flag[8] == '0' and self.bin_flag[11] == '0' else False
        
        #list of mapped bases in reads if read is mapped
        self.mapped_bases = self.find_mapped_elements(self.seq) if self.is_mapped else [] #TO DO: update find_mapped_elements to work without args
        self.mapped_qual = self.find_mapped_elements(self.qual) if self.is_mapped else []

    @property
    def read_range(self) -> tuple[int, int]:
        """
        find the range of mapped read

        Returns:
            returns the start and end of reference positions of mapped read  
        """
        start = self.pos
        end = start + len(self.mapped_bases)

        return start, end
        
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
            clip_idx = cigar_list[-1][1]
            cigar_list = [[l[0], l[1] - clip_idx] for l in cigar_list]
            return cigar_list, clip_idx
        
        return cigar_list, 0

    def find_mapped_elements(self, elements: str) -> list[str]:
        """
        returns a list of mapped elements relative to the reference

        Args:
            elements: sequence or quality scores of mapped read

        Returns:
            list containing only mapped elements
        """

        ele_list = list(elements)
        cigar_list, clip_idx = self.process_cigar()

        #remove soft clipped elements
        if cigar_list[0][0] == 'S':
            del ele_list[:clip_idx]      
        if cigar_list[-1][0] == 'S':
            del ele_list[clip_idx:]

        to_del = []
        #process the element list according to cigar
        for n,l in enumerate(cigar_list):
            if l[0] == 'S':
                continue #already handled
            
            elif l[0] == 'M':
                continue #no action required
            
            elif l[0] == 'D': #add empty string to denote deleted bases - [ele1, ele2, ' ', ele3]
                start_idx = cigar_list[n-1][1]
                stop_idx = l[1]
                diff = stop_idx - start_idx
                ele_list[start_idx:start_idx] = list(' ' * diff)

            elif l[0] == 'I': #combine bases into list to denote insertion - [ele1, [ele2, ele3], ele4]
                start_idx = cigar_list[n-1][1] - 1
                stop_idx = l[1] 
                diff = stop_idx - start_idx - 1
                ele_list[start_idx] = ele_list[start_idx:stop_idx]
                to_del.append([start_idx + 1, diff]) #keep track of index and length to remove bases after combining into list

        # print(f"after processing: {len(ele_list)}, {ele_list}\n")

        #remove elements next to insertions
        for ele in to_del[::-1]: #iterate from end to avoid changing list index
            del ele_list[ele[0]:ele[0]+ele[1]]
        
        return ele_list

    def base_at_pos(self, pos: int) -> str:
        """
        return the base mapped at a position

        Args:
            pos: position relative to the reference

        Returns:
            base at pos
        """
        base = ""
        base_list = list(self.seq)

        # unmapped read
        if not self.is_mapped:
            return ""

        #list conataining mapped bases
        # base_list = self.find_mapped_elements(self.seq)
        base_list = self.mapped_bases

        #check if pos within the range of mapped read
        if self.pos <= pos < self.pos + len(base_list):
            idx = abs(self.pos - pos)
            base = "".join(base_list[idx])
            return "" if base == ' ' else base
        else:
            return ""
    
    def reverse_complement(self, seq: str) -> str:
        comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '':''}
        return "".join([comp_dict[base] for base in seq[::-1]])

    def qual_at_pos(self, pos: int) -> str:
        """
        return the quality score mapped at a position

        Args:
            pos: position relative to the reference

        Returns:
            quality at pos
        """
        qual = ""
        qual_list = list(self.qual)

        #unmapped read
        if not self.is_mapped:
            return ""

        #list conataining quality scores corresponding to mapped bases
        qual_list = self.mapped_qual

        #check if pos within the range of mapped read
        if self.pos <= pos < self.pos + len(qual_list):
            idx = abs(self.pos - pos)
            qual = "".join(qual_list[idx])
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


print(SAM.from_sam("../ex11_data/sams/ERR11767307.sam").primary_reads[9].rname)
sam = SAM.from_sam("../ex11_data/sams/ERR11767307.sam")
for l in sam.reads_at_pos("Fusibacter_paucivorans", 1100):
    print(l.qname, l.seq)
b,q = sam.pileup_at_pos("Fusibacter_paucivorans", 1099)
print(b,q)

# sam_records = [
#     # 'ERR11767307.541398\t163\tFusibacter_paucivorans\t1\t60\t68S83M\t=\t314\t464\tCAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:15\tms:i:136\tAS:i:136\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:105\ts2:i:0\tde:f:0.1807\trl:i:0\n',
#     # 'ERR11767307.1723569\t163\tFusibacter_paucivorans\t1\t60\t60S91M\t=\t376\t529\tAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGT\tF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFF\tNM:i:18\tms:i:146\tAS:i:146\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:117\ts2:i:0\tde:f:0.1978\trl:i:0\n',
#     # 'ERR11767307.67699\t89\tFusibacter_paucivorans\t11\t25\t5S60M26D87M\t=\t11\t0\tCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACATATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGATAACATACCGAAAGGTATGCTAATACGGGATAAT\tFFFFFFFFFFFF-FFFF-FF-FFFFFFFFFF-FFFFFFFF-FFFFFFFFFF-FFFFFFFFFFFFFFFFFFFF-FFFFFFFF-F-FFFFFFFFFFFFFFFFFFFFFFFFFFF55555555555-F55FFFFFFFFFFFFF5F55FFFFFFF-\tNM:i:59\tms:i:222\tAS:i:186\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:45\ts2:i:0\tde:f:0.2237\trl:i:0\n',
#     # 'ERR11767307.11970\t83\tFusibacter_paucivorans\t201\t60\t151M\t=\t1\t-351\tATCTCTTGAATATCAAAGGTGAGCCAGTACAGGATGGACCCGCGTCTGATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCGACGATCAGTAGCCGACCTGAGAGGGTGATCGGCCACATTGGAACTGAGACACGGTCCAAACTCCTAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:20\tms:i:261\tAS:i:263\tnn:i:1\ttp:A:P\tcm:i:11\ts1:i:103\ts2:i:0\tde:f:0.1325\trl:i:0,\n',
#     # 'ERR11767307.1723509\t163\tFusibacter_paucivorans\t107\t56\t5H146M\t=\t514\t558\tACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGATAACATACCGAAAGGTATGCTAATACGGGATAATATATTTGAGAGGCATCTCTTGAATATCAAAGGTGAGCCAGTACAGGATGGACCCGCGTCTGATTAGCTAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFF-FFFFFFFFFFFF-FFFFFFFFFFFFFF-FF-FFFFFFFFFFFFFFFFFF-FFFF\tNM:i:45\tms:i:212\tAS:i:212\tnn:i:0\ttp:A:P\tcm:i:7\ts1:i:68\ts2:i:0\tde:f:0.298\trl:i:0\n',
#     # 'ERR11767307.11970\t163\tFusibacter_paucivorans\t1\t60\t91S60M\t=\t201\t351\tTTAAACAGTAGGTTAATTTATATTAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGC\tFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFF5FFFFFFFF5FFFFFFFFFFFFFFFFF\tNM:i:1\tms:i:118\tAS:i:118\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:103\ts2:i:0\tde:f:0.0167\trl:i:0\n',
#     # 'ERR11767307.1284302\t99\tFusibacter_paucivorans\t1\t60\t23S74M26D54M\t=\t530\t680\tTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:42\tms:i:210\tAS:i:174\tnn:i:0\ttp:A:P\tcm:i:13\ts1:i:101\ts2:i:0\tde:f:0.1318\trl:i:0\n',
#     # 'ERR11767307.12471\t163\tFusibacter_paucivorans\t970\t60\t28M1I14M1I107M\t=\t1397\t556\tCTTGACATCCCAATGACATCTCCTTAATCGGAGAGTTCCCTTCGGGGGCATTGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCTTTAGTTGCCATCATTAAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFF5FFFFFFFFF\tNM:i:27\tms:i:220\tAS:i:220\tnn:i:0\ttp:A:P\tcm:i:12\ts1:i:124\ts2:i:0\tde:f:0.1788\trl:i:0\n',
#     # 'ERR11767307.70912\t99\tFusibacter_paucivorans\t541\t60\t151M\t=\t908\t516\tGGATTTACTGGGCGTAAAGGGTGCGTAGGCGGTCTTTCAAGTCAGGAGTGAAAGGCTACGGCTCAACCGTAGTAAGCTCTTGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFF-FFF\tNM:i:22\tms:i:258\tAS:i:258\tnn:i:0\ttp:A:P\tcm:i:7\ts1:i:100\ts2:i:0\tde:f:0.1457\trl:i:0\n',
#     'ERR11767307.70912\t147\tFusibacter_paucivorans\t908\t60\t101M2I48M\t=\t541\t-516\tGACCCGCACAAGTAGCGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTAAGCTTGACATCCCAATGACATCTCCTTAATCGGAGAGTTCCCTTCGGGGACATTGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTC\tFFFFFFFFF-FFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFF55FFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:30\tms:i:228\tAS:i:226\tnn:i:0\ttp:A:P\tcm:i:10\ts1:i:100\ts2:i:0\tde:f:0.1933\trl:i:0\n'
# ]

# for record in sam_records:
#     read = SAMRead(record)

# pos = 1006
# print(read.qname, pos, read.cigar, read.seq)
# print(repr(read.base_at_pos(pos)))
# print(repr(read.qual_at_pos(pos)))
# print(read.mapped_seq())
# print(read.read_range)
