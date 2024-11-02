#!/usr/bin/env python3

class SAMRead:
    def __init__(self, record:str):
        """
        initialize the fields in record

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

    def base_at_pos():
        pass

    def qual_at_pos():
        pass

    def mapped_seq():
        pass


sam_records = [
    'ERR11767307.541398\t163\tFusibacter_paucivorans\t1\t60\t68S83M\t=\t314\t464\tCAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:15\tms:i:136\tAS:i:136\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:105\ts2:i:0\tde:f:0.1807\trl:i:0\n',
    'ERR11767307.1723569\t163\tFusibacter_paucivorans\t1\t60\t60S91M\t=\t376\t529\tAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGT\tF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFF\tNM:i:18\tms:i:146\tAS:i:146\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:117\ts2:i:0\tde:f:0.1978\trl:i:0\n'
]

for record in sam_records:
    read = SAMRead(record)

print(read.qname, read.flag, (read.bin_flag))
print(read.bin_flag[2], read.bin_flag[4], read.bin_flag[4], read.bin_flag[8])
print(read.is_mapped, read.is_forward, read.is_reverse, read.is_primary)
for n,base in enumerate(read.bin_flag):
    print(n, base)