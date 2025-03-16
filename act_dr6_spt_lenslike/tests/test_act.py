import unittest
import act_dr6_spt_lenslike as apslike
import numpy as np
import os
file_dir = os.path.abspath(os.path.dirname(__file__))
version = apslike.default_version
data_dir = f"{file_dir}/../data/{version}/"


class ACTLikeTest(unittest.TestCase):

    def generic_call(self,variant,lens_only,exp_chisq=None,return_theory=False):
        try:
            ell, cl_tt, cl_ee, cl_bb, cl_te = np.loadtxt(data_dir+'like_corrs/cosmo2017_10K_acc3_lensedCls.dat', unpack=True)
            ellp, _, _, _, _, cl_pp, _, _= np.loadtxt(data_dir+'like_corrs/cosmo2017_10K_acc3_lenspotentialCls.dat', unpack=True)
        except OSError:
            print('Required data file not found at {}'.format(data_file))
            print('Please obtain it and place it correctly.')
            print('The script get-act-data.sh will download and place it.')
            raise

        prefac = 2*np.pi/ell/(ell+1.)
        cl_kk=cl_pp/4*2*np.pi
        cl_bb = cl_bb*prefac
        cl_tt = cl_tt*prefac
        cl_ee = cl_ee*prefac
        cl_te = cl_te*prefac
        data_dict = apslike.load_data(variant,lens_only=lens_only,like_corrections=not(lens_only),version=version)
        ell_kk = ellp
        ell_cmb=ell

        if return_theory:
            chisq,bclkk=apslike.generic_lnlike(data_dict,ell_kk,cl_kk,ell_cmb,cl_tt,cl_ee,cl_te,cl_bb,trim_lmax = 2998,return_theory=True)
            self.assertAlmostEqual(-2*chisq,  exp_chisq, 1)
        else:
            chisq=-2*apslike.generic_lnlike(data_dict,ell_kk,cl_kk,ell_cmb,cl_tt,cl_ee,cl_te,cl_bb,trim_lmax = 2998)
            self.assertAlmostEqual(chisq,  exp_chisq, 1)

    def test_act_baseline_lensonly(self):
        self.generic_call('act_baseline',True,14.06)
    def test_act_baseline(self):
        self.generic_call('act_baseline',False,14.13)
    def test_act_baseline_return_theory(self):
        self.generic_call('act_baseline',False,14.13,return_theory=True)
    def test_actplanck_baseline_lensonly(self):
        self.generic_call('actplanck_baseline',True,21.07)
    def test_actplanck_baseline(self):
        self.generic_call('actplanck_baseline',False,21.46)
    def test_act_extended_lensonly(self):
        self.generic_call('act_extended',True,17.66)
    def test_act_extended(self):
        self.generic_call('act_extended',False,17.72)
    def test_actplanck_extended_lensonly(self):
        self.generic_call('actplanck_extended',True,24.4)
    def test_actplanck_extended(self):
        self.generic_call('actplanck_extended',False,24.75)
    def test_act_polonly_lensonly(self):
        self.generic_call('act_polonly',True,309.71)
    def test_act_cibdeproj_lensonly(self):
        self.generic_call('act_cibdeproj',True,15.16)
    def test_act_cinpaint_lensonly(self):
        self.generic_call('act_cinpaint',True,15.94)
    def test_spt3g_lensonly(self):
        self.generic_call('spt3g',True,19.69)
    # def test_act_spt3g_baseline(self):
    #     self.generic_call('actspt3g_baseline',False,19.69)
    # def test_act_spt3g_extended(self):
    #     self.generic_call('actspt3g_extended',False,26.27)
    def test_act_spt3g_baseline_lensonly(self):
        self.generic_call('actspt3g_baseline',True,33.00)
    def test_act_spt3g_extended_lensonly(self):
        self.generic_call('actspt3g_extended',True,36.53)
    # def test_actplanck_spt3g_baseline(self):
    #     self.generic_call('actplanckspt3g_baseline',False,26.27)
    # def test_actplanck_spt3g_extended(self):
    #     self.generic_call('actplanckspt3g_extended',False,32.85)
    def test_actplanck_spt3g_baseline_lensonly(self):
        self.generic_call('actplanckspt3g_baseline',True,38.17)
    def test_actplanck_spt3g_extended_lensonly(self):
        self.generic_call('actplanckspt3g_extended',True,41.27)

if __name__ == '__main__':
    ACTLikeTest().test_act_baseline_lensonly()
    ACTLikeTest().test_act_baseline()
    ACTLikeTest().test_act_baseline_return_theory()
    ACTLikeTest().test_actplanck_baseline_lensonly()
    ACTLikeTest().test_actplanck_baseline()
    ACTLikeTest().test_act_extended_lensonly()
    ACTLikeTest().test_act_extended()
    ACTLikeTest().test_actplanck_extended_lensonly()
    ACTLikeTest().test_actplanck_extended()
    ACTLikeTest().test_act_polonly_lensonly()
    ACTLikeTest().test_act_cibdeproj_lensonly()
    ACTLikeTest().test_act_cinpaint_lensonly()
    ACTLikeTest().test_spt3g_lensonly()
    # ACTLikeTest().test_act_spt3g_baseline()
    # ACTLikeTest().test_act_spt3g_extended()
    ACTLikeTest().test_act_spt3g_baseline_lensonly()
    ACTLikeTest().test_act_spt3g_extended_lensonly()
    # ACTLikeTest().test_actplanck_spt3g_baseline()
    # ACTLikeTest().test_actplanck_spt3g_extended()
    ACTLikeTest().test_actplanck_spt3g_baseline_lensonly()
    ACTLikeTest().test_actplanck_spt3g_extended_lensonly()