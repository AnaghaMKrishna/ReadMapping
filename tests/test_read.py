import sys
import inspect

import pytest

from .helper_functions import (
    check_class_method,
    check_class_attribute,
    check_attribute_value,
    check_method_output,
)
from . import test_data


try:
    from readmapping.sam_record_processing import SAMRead
except:
    pytest.fail(f"Can't perform remaining tests as imports failed")
    sys.exit()

class TestBasicUsability:
    def test_read_is_a_class(self):
        """Is SAMRead a class"""
        if not inspect.isclass(SAMRead):
            pytest.fail(f"SAMRead is not a class, but an instance of {SAMRead.__class__.__name__}")

    def test_read_class_can_be_instantiated(self):
        """Can the SAMRead class be instantiated from a SAM entry"""
        try:
            _ = SAMRead(test_data.TEST_READ_F_MAPPED)
        except:
            pytest.fail(f"Your SAMRead class __init__ method does not create an instance of the class when given a SAM entry")

class TestReadAttrsExist:
    # Each method in the class gets access to a different instance of this variable to work with using self.read_instance
    read_instance = SAMRead(test_data.TEST_READ_F_MAPPED)
    
    def test_is_mapped_attr_exists(self):
        """Is the 'is_mapped' attribute defined"""
        check_class_attribute(self.read_instance, "is_mapped")

    def test_is_forward_attr_exists(self):
        """Is the 'is_forward' attribute defined"""
        check_class_attribute(self.read_instance, "is_forward")

    def test_is_reverse_attr_exists(self):
        """Is the 'is_reverse' attribute defined"""
        check_class_attribute(self.read_instance, "is_reverse")

    def test_is_primary_attr_exists(self):
        """Is the 'is_primary' attribute defined"""
        check_class_attribute(self.read_instance, "is_primary")

class TestReadMethodsImplemented:
    read_instance = SAMRead(test_data.TEST_READ_F_MAPPED)

    def test_init_method_implemented(self):
        """Is there a callable __init__ method defined"""
        check_class_method(SAMRead, "__init__")

    def test_base_at_pos_method_implemented(self):
        """Is there a callable base_at_pos method defined"""
        check_class_method(SAMRead, "base_at_pos")
    
    def test_qual_at_pos_method_implemented(self):
        """Is there a callable qual_at_pos method defined"""
        check_class_method(SAMRead, "qual_at_pos")

    def test_mapped_seq_method_implemented(self):
        """Is there a callable mapped_seq method defined"""
        check_class_method(SAMRead, "mapped_seq")

class TestReadAttrsBehavior:
    f_mapped_read_instance = SAMRead(test_data.TEST_READ_F_MAPPED)
    r_mapped_read_instance = SAMRead(test_data.TEST_READ_R_MAPPED)
    unmapped_read_instance = SAMRead(test_data.TEST_READ_UNMAPPED)
    secondary_read_instance = SAMRead(test_data.TEST_READ_F_SECONDARY)
    supplemental_read_instance = SAMRead(test_data.TEST_READ_F_SUPPLEMENTAL)

    def test_is_forward_with_forward(self):
        """Is the is_forward attribute set properly"""
        check_attribute_value(
            instance=self.f_mapped_read_instance,
            attr="is_forward",
            expected=True,
            tested_data="forward-mapped read"
        )

    def test_is_reverse_with_forward(self):
        """Is the is_reverse attribute set properly"""
        check_attribute_value(
            instance=self.f_mapped_read_instance,
            attr="is_reverse",
            expected=False,
            tested_data="forward-mapped read"
        )
    
    def test_is_forward_with_reverse(self):
        """Is the is_forward attribute set properly"""
        check_attribute_value(
            instance=self.r_mapped_read_instance,
            attr="is_forward",
            expected=False,
            tested_data="reverse-mapped read"
        )

    def test_is_reverse_with_reverse(self):
        """Is the is_reverse attribute set properly"""
        check_attribute_value(
            instance=self.r_mapped_read_instance,
            attr="is_reverse",
            expected=True,
            tested_data="reverse-mapped read"
        )

    def test_is_mapped_with_mapped(self):
        """Is the is_mapped attribute set properly"""
        check_attribute_value(
            instance=self.f_mapped_read_instance,
            attr="is_mapped",
            expected=True,
            tested_data="mapped read"
        )

    def test_is_mapped_with_unmapped(self):
        """Is the is_mapped attribute set properly"""
        check_attribute_value(
            instance=self.unmapped_read_instance,
            attr="is_mapped",
            expected=False,
            tested_data="unmapped read"
        )

    def test_is_primary_with_primary(self):
        """Is the is_primary attribute set properly"""
        check_attribute_value(
            instance=self.f_mapped_read_instance,
            attr="is_primary",
            expected=True,
            tested_data="primary read mapping"
        )

    def test_is_primary_with_secondary(self):
        """Is the is_primary attribute set properly"""
        check_attribute_value(
            instance=self.secondary_read_instance,
            attr="is_primary",
            expected=False,
            tested_data="secondary read mapping"
        )

    def test_is_primary_with_supplemental(self):
        """Is the is_primary attribute set properly"""
        check_attribute_value(
            instance=self.supplemental_read_instance,
            attr="is_primary",
            expected=False,
            tested_data="supplemental read mapping"
        )

class TestReadAtPosBehavior:
    mapped_read_instance = SAMRead(test_data.TEST_READ_F_MAPPED)
    unmapped_read_instance = SAMRead(test_data.TEST_READ_UNMAPPED)
    non_1_pos_read_instance = SAMRead(test_data.TEST_READ_CIGAR_M_ONLY_NON_1_POS)
    deletion_read_instance = SAMRead(test_data.TEST_READ_CIGAR_D_12)
    insertion_read_instance = SAMRead(test_data.TEST_READ_CIGAR_I_10)
    soft_clip_read_instance = SAMRead(test_data.TEST_READ_SOFT_CLIP)
    hard_clip_read_instance = SAMRead(test_data.TEST_READ_HARD_CLIP)

    def test_init_creates_instance(self):
        """Does the SAMRead.__init__ return an instance of SAMRead"""
        if not isinstance(self.mapped_read_instance, SAMRead):
            pytest.fail(f"Your SAMRead class __init__ method does not create an instance of the class when given a SAM entry")
    
    def test_just_M_read_at_pos(self):
        """Does the base_at_pos method return the right base when the read doesn't map to the first reference base"""
        check_method_output(
            instance=self.non_1_pos_read_instance,
            method="base_at_pos",
            args=(13,), # comma specifies this is a single element tuple. Just here to satisfy syntax rules
            expected="G",
            tested_data="read with only M in CIGAR, but that doesn't map to the start of the reference"
        )
    
    def test_del_read_at_pos(self):
        """Does the base_at_pos method return the right base when there is a deletion in the read at the requested location"""
        check_method_output(
            instance=self.deletion_read_instance,
            method="base_at_pos",
            args=(13,),
            expected="",
            tested_data="read with deletion at requested location"
        )
    
    def test_upstream_del_read_at_pos(self):
        """Does the base_at_pos method return the right base when there is a deletion in the read before the requested location"""
        check_method_output(
            instance=self.deletion_read_instance,
            method="base_at_pos",
            args=(14,),
            expected="G",
            tested_data="read with deletion before requested location"
        )

    def test_ins_read_at_pos(self):
        """Does the base_at_pos method return the right base when there is an insertion in the read at the requested location"""
        check_method_output(
            instance=self.insertion_read_instance,
            method="base_at_pos",
            args=(10,),
            expected="GC",
            tested_data="read with insertion at requested location"
        )

    def test_upstream_ins_read_at_pos(self):
        """Does the base_at_pos method return the right base when there is an insertion in the read before the requested location"""
        check_method_output(
            instance=self.insertion_read_instance,
            method="base_at_pos",
            args=(14,),
            expected="G",
            tested_data="read with insertion before (i.e., 5') requested location"
        )
    
    def test_downstream_del_read_at_pos(self):
        """Is the base_at_pos method impacted if there is a deletion in the read after the requested location"""
        check_method_output(
            instance=self.deletion_read_instance,
            method="base_at_pos",
            args=(10,),
            expected="G",
            tested_data="read with deletion after (i.e., 3') requested location"
        )

    def test_mapped_right_of_requested_pos_read_at_pos(self):
        """Does base_at_pos return an empty string for reads that map to the right of requested location"""
        check_method_output(
            instance=self.mapped_read_instance,
            method="base_at_pos",
            args=(10,),
            expected="",
            tested_data="read maps to the right of requested requested location"
        )

    def test_first_base_in_read_read_at_pos(self):
        """Does base_at_pos return an empty string for reads that map to the right of requested location"""
        check_method_output(
            instance=self.mapped_read_instance,
            method="base_at_pos",
            args=(16,),
            expected="G",
            tested_data="requested is the first base in the read"
        )
    
    def test_mapped_left_of_requested_pos_read_at_pos(self):
        """Does base_at_pos return an empty string for reads that map to the left of requested location"""
        check_method_output(
            instance=self.mapped_read_instance,
            method="base_at_pos",
            args=(167,),
            expected="",
            tested_data="read maps to the left of requested requested location"
        )

    def test_last_base_in_read_read_at_pos(self):
        """Does base_at_pos return the right read when the requested pos is the last base in the read"""
        check_method_output(
            instance=self.mapped_read_instance,
            method="base_at_pos",
            args=(166,),
            expected="C",
            tested_data="requested is the last base in the read"
        )

    def test_unmapped_read_pos_read_at_pos(self):
        """Does the base_at_pos method return an empty string for unmapped reads"""
        check_method_output(
            instance=self.unmapped_read_instance,
            method="base_at_pos",
            args=(10,),
            expected="",
            tested_data="unmapped read"
        )

    def test_soft_clipped_upstream_read_at_pos(self):
        """Does the base_at_pos method correctly handle soft clipping"""
        check_method_output(
            instance=self.soft_clip_read_instance,
            method="base_at_pos",
            args=(13,),
            expected="G",
            tested_data="read with soft clipping before requested location"
        )
    
    def test_hard_clipped_upstream_read_at_pos(self):
        """Does the base_at_pos method correctly handle hard clipping"""
        check_method_output(
            instance=self.hard_clip_read_instance,
            method="base_at_pos",
            args=(13,),
            expected="G",
            tested_data="read with hard clipping before requested location"
        )

    def test_internally_mapped_read_at_pos(self):
        """Does the base_at_pos method return the right base for a position in the middle of the reference"""
        check_method_output(
            instance=self.non_1_pos_read_instance,
            method="base_at_pos",
            args=(50,),
            expected="C",
            tested_data="read with only M in CIGAR, that maps to somewhere in the middle of the reference"
        )

class TestMappedSeqBehavior:
    middle_mapped_read_instance = SAMRead(test_data.TEST_READ_CIGAR_M_ONLY_NON_1_POS)
    deletion_read_instance = SAMRead(test_data.TEST_READ_CIGAR_D_12)
    insertion_read_instance = SAMRead(test_data.TEST_READ_CIGAR_I_10)
    soft_clipped_read_instance = SAMRead(test_data.TEST_READ_SOFT_CLIP)
    hard_clipped_read_instance = SAMRead(test_data.TEST_READ_HARD_CLIP)
    unmapped_read_instance = SAMRead(test_data.TEST_READ_UNMAPPED)

    def test_internally_mapped_mapped_seq(self):
        """Does the mapped_seq method return the right sequence for a read mapped to the middle of the reference"""
        check_method_output(
            instance=self.middle_mapped_read_instance,
            method="mapped_seq",
            expected="CAGCAGCCGCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTGAGCTAGAGTACGGTAGAGGGT",
            tested_data="read with only M in CIGAR, that maps to somewhere in the middle of the reference"
        )

    def test_deletion_mapped_seq(self):
        """Does the mapped_seq method handle deletions correctly"""
        check_method_output(
            instance=self.deletion_read_instance,
            method="mapped_seq",
            expected="GTGCCAGCAGCC-GCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTNGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTGAGCTAGAGTACGGTAG",
            tested_data="read with a deletion relative to the reference"
        )
    
    def test_insertion_mapped_seq(self):
        """Does the mapped_seq method handle insertions correctly"""
        check_method_output(
            instance=self.insertion_read_instance,
            method="mapped_seq",
            expected="GTGCCAGCAGCCGCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTNGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTGAGCTAGAGTACGGTAG",
            tested_data="read with an insertion relative to the reference"
        )

    def test_soft_clipped_mapped_seq(self):
        """Does the mapped_seq method handle soft clipping correctly"""
        check_method_output(
            instance=self.soft_clipped_read_instance,
            method="mapped_seq",
            expected="GTGCCAGCAGCCGCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTG",
            tested_data="read with an insertion relative to the reference"
        )

    def test_hard_clipped_mapped_seq(self):
        """Does the mapped_seq method handle hard clipping correctly"""
        check_method_output(
            instance=self.hard_clipped_read_instance,
            method="mapped_seq",
            expected="GTGCCAGCAGCCGCGGTAATACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTG",
            tested_data="read with an insertion relative to the reference"
        )
