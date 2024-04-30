import unittest
import pandas as pd
import numpy as np

from casedi.survival import ase_compare_survival

example_ase_table = pd.DataFrame({'gene': np.repeat('A',20),
                                 'ase': np.concatenate([np.repeat(0,10),np.repeat(1,10)]),
                                 'sampleID': np.arange(1,21,1) # mismatching sampleID
                                 })
example_survival_table = pd.DataFrame({'survival_time': np.concatenate([np.repeat(10,10),np.repeat(2,10)]),
                                       'survival_status': np.concatenate([np.repeat(1,10),np.repeat(0,10)]),
                                      })
example_expression_table = pd.DataFrame({'A': np.arange(1,21,1)})

example_expression_table2 = pd.DataFrame({'A': np.arange(1,22,1)})

class TestSurvival(unittest.TestCase):
    """Test `casedi.survival` functions"""
    
    def test_mismatch_id(self):
        with self.assertRaises(Exception) as context:
            ase_compare_survival(gene='A',
                                 ase_table=example_ase_table,
                               survival_table=example_survival_table,
                               expression_table=example_expression_table,
                               stratify_expression='ase high')

        self.assertTrue('Sample IDs of tables do not match' in str(context.exception))
        
    def test_no_expr(self):
        with self.assertRaises(ValueError) as context:
            ase_compare_survival(gene='A',
                                 ase_table=example_ase_table,
                               survival_table=example_survival_table,
                               stratify_expression='ase high')

        self.assertTrue("Expression table must be passed if stratify_expression is set" in str(context.exception))
    
    def test_mismatch_length(self):
        with self.assertRaises(Exception) as context:
            ase_compare_survival(gene='A',
                                 ase_table=example_ase_table,
                               survival_table=example_survival_table,
                               expression_table=example_expression_table2,
                               stratify_expression='ase high')

        self.assertTrue("Unequal number of Sample IDs in tables" in str(context.exception))

if __name__ == '__main__':
    unittest.main()