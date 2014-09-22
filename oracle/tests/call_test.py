#dodgy test script

def kek():
    import os
    os.system('bench.py')

import unittest

class TestUM(unittest.TestCase):
 
    def setUp(self):
        pass

    def test_strings_a_3(self):
        self.assertIsNone(kek())
        print kek()
        
if __name__ == '__main__':
    unittest.main()
