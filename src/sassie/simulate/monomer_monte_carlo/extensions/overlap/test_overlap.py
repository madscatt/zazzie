import sys
import unittest
import numpy as np
sys.path.append('./build/lib.linux-x86_64-3.9/')
import overlap

class TestOverlap(unittest.TestCase):
    def test_overlap_far(self):
        # Create a 2D numpy array with 3 columns
        array = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
        cut = 1.75

        # Call the overlap function
        result = overlap.overlap(array, cut)
        #   5.196152422706632

        # Add your assertions here
        # For example, if overlap is supposed to return True if any pair of points is closer than `cut`:
        self.assertEqual(result, False)

    def test_overlap_close(self):
        # Create a 2D numpy array with 3 columns
        # Points are closer than 1.5
        array = np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]], dtype=np.float32)
        cut = 1.75

        # Call the overlap function
        result = overlap.overlap(array, cut)
        # 1.7320508075688772

        # Add your assertions here
        # For example, if overlap is supposed to return True if any pair of points is closer than `cut`:
        self.assertEqual(result, True)

if __name__ == '__main__':
    unittest.main()
