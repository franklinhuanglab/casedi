import unittest
import pandas as pd
import numpy as np

from casedi_test.expression import ase_compare_expression

example_ase_table = pd.DataFrame({'gene': np.repeat('A',20),
                                 'ase': np.concatenate([np.repeat(0,10),np.repeat(1,10)]),
                                 'sampleID': np.arange(0,20,1) 
                                 })
example_ase_table2 = pd.DataFrame({'gene': np.repeat('A',20),
                                 'ase': np.repeat(0,20), # no ASE samples
                                 'sampleID': np.arange(0,20,1) 
                                 })
example_expression_table = pd.DataFrame({'A': np.arange(0,20,1)}) 
example_expression_table.index = np.arange(1,21,1) # mismatching sampleID

example_expression_table2 = pd.DataFrame({'A': np.arange(1,22,1)}) # mismatching length

example_expression_table3 = pd.DataFrame({'A': np.arange(1,21,1)}) 


class TestExpression(unittest.TestCase):
    """Test `casedi.expresion` functions"""
    def test_mismatch_length(self):
        with self.assertRaises(Exception) as context:
            ase_compare_expression(gene='A',
                                   ase_table=example_ase_table,
                                   expression_table=example_expression_table2,
                                   )

        self.assertTrue("Unequal number of Sample IDs in tables" in str(context.exception))

    def test_mismatch_id(self):
        with self.assertRaises(Exception) as context:
            ase_compare_expression(gene='A',
                                   ase_table=example_ase_table,
                                   expression_table=example_expression_table,
                                   )

        self.assertTrue('Sample IDs of tables do not match' in str(context.exception))
        
    def test_no_ase(self):
        with self.assertRaises(Exception) as context:
            ase_compare_expression(gene='A',
                                   ase_table=example_ase_table2,
                                   expression_table=example_expression_table3)

        self.assertTrue("No samples showing ASE" in str(context.exception))
    


if __name__ == '__main__':
    unittest.main()