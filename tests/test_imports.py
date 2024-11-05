import pytest

# Grouping related tests into a class is just for organization purposes
class TestImports:
    def test_magnumopus_import(self):
        """Can the readmapping package be imported"""
        try:
            import readmapping
        except:
            pytest.fail(f"Your readmapping package could not be imported")
    
    def test_read_module_import(self):
        """Can the read module be imported"""
        try:
            from readmapping import sam_record_processing
        except:
            pytest.fail(f"Your sam_record_processing module could not be imported")
    
    def test_read_class_import(self):
        """Can the SAMRead class be imported"""
        try:
            from readmapping.sam_record_processing import SAMRead
        except:
            pytest.fail(f"Your SAMRead class could not be imported")
