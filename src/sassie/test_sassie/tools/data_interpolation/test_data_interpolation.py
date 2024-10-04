import os
import shutil
import filecmp
from unittest import main, TestCase
from unittest.mock import patch

import sassie.tools.data_interpolation.gui_mimic_data_interpolation as gui_mimic_data_interpolation

pdb_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'pdb_common') + os.path.sep
dcd_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'dcd_common') + os.path.sep
other_data_path = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', '..', 'data', 'other_common') + os.path.sep
module_data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'tools', 'data_interpolation') + os.path.sep

paths = {'pdb_data_path': pdb_data_path, 'dcd_data_path': dcd_data_path,
         'other_data_path': other_data_path, 'module_data_path': module_data_path}


class Test_Data_Interpolation(TestCase):

    module = 'data_interpolation'

    def setUp(self):
        gui_mimic_data_interpolation.test_variables(self, paths)

    @patch('sassie.tools.data_interpolation.gui_mimic_data_interpolation.run_module')
    def test_1(self, MockClass):
        '''
        test SANS data file
        '''
        instance = MockClass.return_value
        instance.some_method.return_value = 'expected_value'

        result = gui_mimic_data_interpolation.run_module(self)

        # confirm output data file is correct
        outfile = os.path.join(self.run_name, self.module, self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'sans_data.dat')
        self.assertEqual(filecmp.cmp(outfile, correct_outfile), True)

        # confirm that output stn_data file is correct
        outfile = os.path.join(self.run_name, self.module,
                               'stn_'+self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'stn_sans_data.dat')
        self.assertEqual(filecmp.cmp(outfile, correct_outfile), True)

    '''
    def test_2(self):
        """
        test SAXS data file
        """
        self.expdata = os.path.join(other_data_path, 'trunc2a_saxs.sub')
        self.ofile = 'trunc2a.dat'
        self.io = '0.031'
        self.dq = '0.007'
        self.maxpoints = '72'

        gui_mimic_data_interpolation.run_module(self)

        # confirm output data file is correct
        outfile = os.path.join(self.run_name, self.module, self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'trunc2a.dat')
        self.assertEqual(filecmp.cmp(outfile, correct_outfile), True)

        # confirm that output stn_data file is correct
        outfile = os.path.join(self.run_name, self.module,
                               'stn_'+self.ofile)
        correct_outfile = os.path.join(
            module_data_path, 'stn_trunc2a.dat')
        self.assertEqual(filecmp.cmp(outfile, correct_outfile), True)
    '''

    def tearDown(self):
        if os.path.exists(self.run_name):
            shutil.rmtree(self.run_name)


if __name__ == '__main__':
    main()