#!/usr/bin/env python3

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

        self.bin_flag = f"{self.flag:012b}"[::-1] #reverse to get LSB at index 0
    
        self.is_mapped = True if self.bin_flag[2] == '0' else False
        self.is_forward = True if self.bin_flag[4] == '0' else False
        self.is_reverse = False if self.bin_flag[4] == '0' else True
        self.is_primary = True if self.bin_flag[8] == '0' and self.bin_flag[11] == '0' else False

    def base_at_pos(self, pos: int) -> str:
        """
        return the base mapped at a position

        Args:
            pos: position relative to the reference

        Returns:
            str: base at specific position
        """
        num = ''
        prev_idx = 0
        cigar_list = []
        #read in CIGAR string to get the read index and reference index - [character, read_index, reference_index]
        for c in self.cigar:
            if c not in "MIDSH":
                num += c
            else:
                cigar_list.append([c, prev_idx + int(num), (self.pos - 1) + prev_idx + int(num)])
                prev_idx = cigar_list[-1][1]
                num = ''

        print(cigar_list) 

        base = ''
        soft_offset = 0
        #process the cigar_list to offset index for soft clipped regions
        if cigar_list[0][0] == 'S':
            soft_offset = cigar_list[0][1]
            for l in cigar_list:
                l[1] -= soft_offset
                l[2] -= soft_offset

        print(cigar_list)

        base = ''
        for l in cigar_list:
            if self.pos <= pos <= l[2]: #check if pos is within the read range
                if l[0] == 'M':
                    seq_idx = abs(pos + soft_offset - self.pos)
                    print(f"seq_idx: {seq_idx}")
                    base = self.seq[seq_idx] if self.is_forward else self.reverse_complement(self.seq[seq_idx])
                    break
                elif l[0] == 'D':
                    break
                elif l[0] == 'I':
                    pass

        return base

    def reverse_complement(self, seq: str) -> str:
        comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '':''}
        return "".join([comp_dict[base] for base in seq[::-1]])

    def qual_at_pos():
        pass

    def mapped_seq():
        pass


sam_records = [
    # 'ERR11767307.541398\t163\tFusibacter_paucivorans\t1\t60\t68S83M\t=\t314\t464\tCAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:15\tms:i:136\tAS:i:136\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:105\ts2:i:0\tde:f:0.1807\trl:i:0\n',
    # 'ERR11767307.1723569\t163\tFusibacter_paucivorans\t1\t60\t60S91M\t=\t376\t529\tAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGT\tF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFF\tNM:i:18\tms:i:146\tAS:i:146\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:117\ts2:i:0\tde:f:0.1978\trl:i:0\n',
    # 'ERR11767307.67699\t89\tFusibacter_paucivorans\t11\t25\t5S60M26D87M\t=\t11\t0\tCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACATATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGATAACATACCGAAAGGTATGCTAATACGGGATAAT\tFFFFFFFFFFFF-FFFF-FF-FFFFFFFFFF-FFFFFFFF-FFFFFFFFFF-FFFFFFFFFFFFFFFFFFFF-FFFFFFFF-F-FFFFFFFFFFFFFFFFFFFFFFFFFFF55555555555-F55FFFFFFFFFFFFF5F55FFFFFFF-\tNM:i:59\tms:i:222\tAS:i:186\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:45\ts2:i:0\tde:f:0.2237\trl:i:0\n',
    # 'ERR11767307.11970\t83\tFusibacter_paucivorans\t201\t60\t151M\t=\t1\t-351\tATCTCTTGAATATCAAAGGTGAGCCAGTACAGGATGGACCCGCGTCTGATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCGACGATCAGTAGCCGACCTGAGAGGGTGATCGGCCACATTGGAACTGAGACACGGTCCAAACTCCTAC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:20\tms:i:261\tAS:i:263\tnn:i:1\ttp:A:P\tcm:i:11\ts1:i:103\ts2:i:0\tde:f:0.1325\trl:i:0,\n',
    # 'ERR11767307.1723509\t163\tFusibacter_paucivorans\t107\t56\t151M\t=\t514\t558\tGGCGGACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGATAACATACCGAAAGGTATGCTAATACGGGATAATATATTTGAGAGGCATCTCTTGAATATCAAAGGTGAGCCAGTACAGGATGGACCCGCGTCTGATTAGCTAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFF-FFFFFFFFFFFF-FFFFFFFFFFFFFF-FF-FFFFFFFFFFFFFFFFFF-FFFF\tNM:i:45\tms:i:212\tAS:i:212\tnn:i:0\ttp:A:P\tcm:i:7\ts1:i:68\ts2:i:0\tde:f:0.298\trl:i:0\n',
    # 'ERR11767307.11970\t163\tFusibacter_paucivorans\t1\t60\t91S60M\t=\t201\t351\tTTAAACAGTAGGTTAATTTATATTAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGC\tFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFF5FFFFFFFF5FFFFFFFFFFFFFFFFF\tNM:i:1\tms:i:118\tAS:i:118\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:103\ts2:i:0\tde:f:0.0167\trl:i:0\n',
    # 'ERR11767307.1284302\t99\tFusibacter_paucivorans\t1\t60\t23S74M26D54M\t=\t530\t680\tTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGTGAGTAACGCGTGGGTAACCTACCCTGTACACACGGAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:42\tms:i:210\tAS:i:174\tnn:i:0\ttp:A:P\tcm:i:13\ts1:i:101\ts2:i:0\tde:f:0.1318\trl:i:0\n',
    'ERR11767307.12471\t163\tFusibacter_paucivorans\t970\t60\t28M1I14M1I107M\t=\t1397\t556\tCTTGACATCCCAATGACATCTCCTTAATCGGAGAGTTCCCTTCGGGGGCATTGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCTTTAGTTGCCATCATTAAG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFFFFFFFFFF5FFFFFFFFF\tNM:i:27\tms:i:220\tAS:i:220\tnn:i:0\ttp:A:P\tcm:i:12\ts1:i:124\ts2:i:0\tde:f:0.1788\trl:i:0\n'
]

for record in sam_records:
    read = SAMRead(record)

print(read.qname, read.cigar, read.seq)
# print(read.bin_flag[2], read.bin_flag[4], read.bin_flag[4], read.bin_flag[8])
# print(read.is_mapped, read.is_forward, read.is_reverse, read.is_primary)
print(read.base_at_pos(997))
