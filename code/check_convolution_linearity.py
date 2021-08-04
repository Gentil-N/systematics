"""IMPORTS"""

import logging as log
log.basicConfig(level=log.INFO)

from core import SampleTime

from math import cos
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math
import cmath
import toast
import toast.todmap
import toast.pipeline_tools
from toast.mpi import MPI
from toast_litebird.hardware import Hardware
from toast_litebird.focalplane import sim_telescope_detectors
from litebird_sim import Imo
import litebird_sim as lbs
log.info("python import succesfull")

"""INPUTS"""

NSIDE = 64
NSIDE_HIGH = NSIDE
LMAX_HIGH = NSIDE_HIGH * 2
ELL_MAX = 10
SAMPLE_TIME = SampleTime(0, 0, 1, 0, 0) # One hour of convolution
EDGE = 10e-16
OUTDIR = "./output/"

"""CONSTANTS"""

NUM_PIX = hp.nside2npix(NSIDE)
MAP_ZERO = np.zeros(NUM_PIX)
ARRAY_ZERO_ALM = hp.map2alm(MAP_ZERO, lmax = ELL_MAX)

"""FUNCS"""

def get_complex_from_ang(theta):
       return complex(math.cos(math.radians(theta)), math.sin(math.radians(theta)))

def create_alm_coeff(ell: int, emm: int, mode: float):
       index = hp.Alm.getidx(lmax = ELL_MAX, l = ell, m = emm)
       alm_coeff = ARRAY_ZERO_ALM.copy()
       alm_coeff[index] = get_complex_from_ang(mode)
       return alm_coeff

def create_map_from_alm_coeff(alm_coeff: np.ndarray):
       return hp.alm2map(alm_coeff, nside = NSIDE, lmax = ELL_MAX)

def create_alm_coeff_with_map(ell: int, emm: int, mode: float):
       alm_coeff = create_alm_coeff(ell, emm, mode)
       return (create_map_from_alm_coeff(alm_coeff), alm_coeff)

def write_blm_from_map(in_map: np.ndarray, out_file: str):
       bl, blm = hp.anafast(np.array([in_map,in_map,in_map]), lmax=LMAX_HIGH, iter=0, alm=True)
       hp.write_alm(out_file, blm, overwrite=True)

def check_map_differences(map_a: np.ndarray, map_b: np.ndarray):
       array_size = min(map_a.size, map_b.size)
       for i in range(array_size):
              diff = abs(map_a[i] - map_b[i])
              if diff > EDGE:
                     return False
       return True

def convolve(comm: toast.Comm, data: toast.Data, nside: int, nside_high: int, beam_file: str, sky_file: str, outdir: str, outprefix: str):

       toast.todmap.OpPointingHpix(nside=nside, nest=True, mode="IQU").exec(data)

       name = "signal"
       toast.tod.OpCacheClear(name).exec(data)
       conviqt = toast.todmap.OpSimConviqt(
           comm.comm_rank,
           sky_file,
           beam_file,
           lmax=2*nside_high,
           beammmax=2*nside_high,
           pol=True,
           fwhm=0,
           order=13,
           calibrate=True,
           dxx=False,
           out=name,
           quat_name=None,
           flag_name=None,
           flag_mask=255,
           common_flag_name=None,
           common_flag_mask=255,
           apply_flags=False,
           remove_monopole=False,
           remove_dipole=False,
           normalize_beam=True,
           verbosity=1,
       )
       conviqt.exec(data)

       mapmaker = toast.todmap.OpMapMaker(
           nside=nside,
           nnz=3,
           outdir=outdir,
           outprefix=outprefix,
           baseline_length=10,
           iter_max=100,
           use_noise_prior=False,
       )
       mapmaker.exec(data)

       binned_map = hp.read_map(outdir + outprefix + "binned.fits") # Get the 'binned' map
       return binned_map

def create_blm_files(a: float, ell_i: int, emm_i: int, mode_i: float, b: float, ell_j: int, emm_j: int, mode_j: float, map_i_name: str, map_j_name: str, total_map_name: str):
       map_i = create_alm_coeff_with_map(ell_i, emm_i, mode_i)[0]
       map_j = create_alm_coeff_with_map(ell_j, emm_j, mode_j)[0]
       total_map = a * map_i + b * map_j
       write_blm_from_map(map_i, OUTDIR + map_i_name)
       write_blm_from_map(map_j, OUTDIR + map_j_name)
       write_blm_from_map(total_map, OUTDIR + total_map_name)

"""MAIN"""

# MPI Loads
mpiworld, procs, rank = toast.mpi.get_world()
comm = toast.mpi.Comm(mpiworld)
# Get scanning strategy
imo = Imo(flatfile_location="./data")
scanning_strategy = lbs.SpinningScanningStrategy.from_imo(imo, "5e3bf49f-53f3-4f80-a65c-6191a4820b86")
# Get detectors
hw = Hardware("./data/hardware_config_files/20Hardware.toml.gz")
log.info("raw file data loaded")

class args:
    sample_rate = 19 # Frame rate (per second)
    hwp_rpm = None
    hwp_step_deg = None
    hwp_step_time_s = None
    spin_period_min = 60*scanning_strategy.spin_rate_hz
    spin_angle_deg = 50
    prec_period_min = 60*scanning_strategy.precession_rate_hz
    prec_angle_deg = 45 
    coord = "E"
    nside = 64
    nnz = 3
# Focalplane = dictionary : [detector (key)] -> [quaternion tuple (value)]
focalplane = {list(hw.data['detectors'].items())[0][0]:list(hw.data['detectors'].items())[0][1]}
# Get detector list
detectors = sorted(focalplane.keys())
# Get quaterion tuple list
detquats = {}
for d in detectors:
    detquats[d] = focalplane[d]["quat"]

nsample = int(SAMPLE_TIME.get_total_seconds() * args.sample_rate)
log.info("detector's infos created")

start_sample = 0
start_time = 0
iobs = 0

tod = toast.todmap.TODSatellite(
        comm.comm_group,
        detquats,
        nsample,
        firstsamp=start_sample,
        firsttime=start_time,
        rate=args.sample_rate,
        spinperiod=args.spin_period_min,
        spinangle=args.spin_angle_deg,
        precperiod=args.prec_period_min,
        precangle=args.prec_angle_deg,
        coord=args.coord,
        #detranks=comm.group_size,
        hwprpm=args.hwp_rpm,
        hwpstep=args.hwp_step_deg,
        hwpsteptime=args.hwp_step_time_s,
    )

# Create empty array with all data
precquat = np.empty(4 * tod.local_samples[1], dtype=np.float64).reshape((-1, 4))
# Constantly slewing precession axis
toast.todmap.slew_precession_axis(
        precquat, # Result
        firstsamp=start_sample + tod.local_samples[0], # Begining sample
        samplerate=args.sample_rate, # Frame rate (per second)
        degday=360.0 / 365.25 # Rotation rate (deg per day)
        )
tod.set_prec_axis(qprec=precquat)

obs = {}
obs["name"] = "science_{:05d}".format(iobs)
obs["tod"] = tod
obs["intervals"] = None
obs["baselines"] = None
obs["id"] = iobs
# Create Focalplane from a dictionary
obs["focalplane"] = toast.pipeline_tools.Focalplane(detector_data=focalplane)
log.info("observation's data created")

data = toast.Data(comm)
data.obs.append(obs) # 'append' is necessray (don't equal)


a = 1.0
b = 1.0
create_blm_files(a, 1, 0, 0.0, b, 1, 1, 0.0, "map_i.fits", "map_j.fits", "total_map.fits")
binned_map_i = convolve(comm, data, NSIDE, NSIDE_HIGH, "map_i.fits", "./data/maptotlm_freq0.fits", OUTDIR, "map_i_")
binned_map_j = convolve(comm, data, NSIDE, NSIDE_HIGH, "map_j.fits", "./data/maptotlm_freq0.fits", OUTDIR, "map_j_")
binned_total_map = convolve(comm, data, NSIDE, NSIDE_HIGH, "total_map.fits", "./data/maptotlm_freq0.fits", OUTDIR, "total_map_")

binned_test_linearity = a * binned_map_i + b * binned_map_j

plt.figure(num = "Results")
hp.mollview(binned_test_linearity, sub=[1, 2, 1], title = "Check")
hp.mollview(binned_test_linearity, sub=[1, 2, 2], title = "Good")
plt.show()

print(check_map_differences(binned_test_linearity, binned_total_map))