#!/usr/bin/env python3

class SAMReads:
    def __init__(self, record):
        pass
    
    is_mapped = None
    is_forward = None
    is_reverse = None
    is_primary = None

    def base_at_pos():
        pass

    def qual_at_pos():
        pass

    def mapped_seq():
        pass


sam_records = [
    'ERR11767307.541398\t163\tFusibacter_paucivorans\t1\t60\t68S83M\t=\t314\t464\tCAAGAAACAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNM:i:15\tms:i:136\tAS:i:136\tnn:i:0\ttp:A:P\tcm:i:8\ts1:i:105\ts2:i:0\tde:f:0.1807\trl:i:0\n'
    'ERR11767307.1723569\t163\tFusibacter_paucivorans\t1\t60\t60S91M\t=\t376\t529\tAAACCATAAAGCCAGATATTTTGATAACAATAGTATCTGAGCCTGATAAACTTTTATTTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTTGAGCGATTTACTTCGGTAAAGAGCGGCGGACGGGT\tF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5FFFFFFFFFFFFF\tNM:i:18\tms:i:146\tAS:i:146\tnn:i:0\ttp:A:P\tcm:i:9\ts1:i:117\ts2:i:0\tde:f:0.1978\trl:i:0\n'
]

for record in sam_records:
    read = SAMReads(record)
