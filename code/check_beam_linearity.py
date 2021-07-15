"""IMPORTS"""

from math import cos
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import cmath

"""INPUTS"""

NSIDE = 64
ELL_MAX = 10
EDGE = 10e-16
ENABLE_PRINT = True
ENABLE_PLOT = False

"""CONSTANTS"""

NUM_PIX = hp.nside2npix(NSIDE)
MAP_ZERO = np.zeros(NUM_PIX)
ARRAY_ZERO_ALM = hp.map2alm(MAP_ZERO, lmax = ELL_MAX)

"""FUNCS"""

def get_complex_from_ang(theta):
       return complex(math.cos(math.radians(theta)), math.sin(math.radians(theta)))

def create_alm_coeff(ell: int, emm: int):
       index = hp.Alm.getidx(lmax = ELL_MAX, l = ell, m = emm)
       alm_coeff = ARRAY_ZERO_ALM.copy()
       alm_coeff[index] = get_complex_from_ang(0.0)
       return alm_coeff

def create_map_from_alm_coeff(alm_coeff: np.ndarray):
       return hp.alm2map(alm_coeff, nside = NSIDE, lmax = ELL_MAX)

def create_alm_coef_with_map(ell: int, emm: int):
       alm_coeff = create_alm_coeff(ell, emm)
       return (create_map_from_alm_coeff(alm_coeff), alm_coeff)

def check_map_differences(map_a: np.ndarray, map_b: np.ndarray):
       array_size = min(map_a.size, map_b.size)
       for i in range(array_size):
              diff = abs(map_a[i] - map_b[i])
              if diff > EDGE:
                     return False
       return True

def check_linearity_between_emm_for_ell(base_ell: int, emm_a: int, emm_b: int):
       alm_coeff_a = create_alm_coeff(base_ell, emm_a)
       alm_coeff_b = create_alm_coeff(base_ell, emm_b)
       alm_coeff_combined = alm_coeff_a + alm_coeff_b
       map_a = create_map_from_alm_coeff(alm_coeff_a)
       map_b = create_map_from_alm_coeff(alm_coeff_b)
       map_combined = create_map_from_alm_coeff(alm_coeff_combined)
       map_check_combined = map_a + map_b
       check = check_map_differences(map_combined, map_check_combined)
       return (emm_a, emm_b, map_a, map_b, map_combined, map_check_combined, check)

def check_linearity_for_ell(ell: int):
       if ENABLE_PRINT:
              print("Current ell : %s"%ell)
       check_count = int((ell + 1) * ell / 2)
       if check_count == 0:
              return
       map_list = []
       for emm_a in range(ell + 1):
              for emm_b in range(emm_a + 1, ell + 1):
                     if ENABLE_PRINT:
                            print("Checking between (%s, %s) : "%(emm_a, emm_b), end = '')
                     map_list.append(check_linearity_between_emm_for_ell(ell, emm_a, emm_b))
                     if ENABLE_PRINT:
                            print(map_list[-1][6])
       if ENABLE_PLOT:
              row = check_count
              col = 4
              plt.figure(num = "Results for ell : %s"%ell)
              for i in range(check_count):
                     hp.mollview(map_list[i][2], sub=[row, col, i * 4 + 1], title="Emm : %s"%map_list[i][0])
                     hp.mollview(map_list[i][3], sub=[row, col, i * 4 + 2], title="Emm : %s"%map_list[i][1])
                     hp.mollview(map_list[i][4], sub=[row, col, i * 4 + 3], title="Good")
                     hp.mollview(map_list[i][5], sub=[row, col, i * 4 + 4], title="Check")
              plt.show()

"""MAIN"""

for i in range(ELL_MAX + 1) :
       check_linearity_for_ell(i)