import random
import sys
import unittest
from pathlib import Path

PROJECT_PATH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_PATH / "src"))

from config_code.complements.shuffled import shuffled

class ShuffleTestCase(unittest.TestCase):
    def test_shuffle(self):
        file_path = PROJECT_PATH / (
            "src/config_code/hpv_step2/targets/"
            "HPV16_E6starI_E7_spliced_region.rna.fasta"
        )
        
        # Verify the file exists before running the test (optional, but good practice)
        self.assertTrue(file_path.exists(), "FASTA file path does not exist.")
        
        sequence = "ACGUACGU"
        result = shuffled(sequence, random.Random(42))
        self.assertTrue(result)
        self.assertCountEqual(result, sequence)

if __name__ == "__main__":
    unittest.main()
