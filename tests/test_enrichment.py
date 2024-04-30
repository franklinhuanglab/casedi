import unittest
import pandas as pd
import numpy as np

from casedi.enrichment import ase_compare_groups

example_ase_table = pd.DataFrame({'gene': np.repeat('A',20),
                                 'ase': np.concatenate([np.repeat(0,10),np.repeat(1,10)]),
                                 'sampleID': np.arange(0,20,1) 
                                 })

group1 = np.array([])
group2 = np.arange(0,10,1) 
group3 = np.arange(13,25,1) 

class TestEnrichment(unittest.TestCase):
    """Test `casedi.enrichment` functions"""
    def test_empty_group(self):
        with self.assertRaises(Exception) as context:
            ase_compare_groups(ase_table=example_ase_table,
                               group1_ids=group1,
                               group2_ids=group2
                              )

        self.assertTrue("One or more groups is empty" in str(context.exception))

    def test_wrong_ids(self):
        with self.assertRaises(Exception) as context:
            ase_compare_groups(ase_table=example_ase_table,
                               group1_ids=group2,
                               group2_ids=group3
                              )

        self.assertTrue('Some group1 or group2 ids were not found in ASE table' in str(context.exception))

if __name__ == '__main__':
    unittest.main()