"""HELPERS"""

class SampleTime:
       seconds: int
       minutes: int
       hours: int
       days: int
       years: int
       def __init__(self, seconds: int = 0, minutes: int = 0, hours: int = 0, days: int = 0, years: int = 0) :
              self.seconds = seconds
              self.minutes = minutes
              self.hours = hours
              self.days = days
              self.years = years

def get_total_sample_time_seconds(sample_time: SampleTime) :
       seconds_per_minutes = 60
       seconds_per_hours = 60 * seconds_per_minutes
       seconds_per_day = 24 * seconds_per_hours
       seconds_per_year = 365 * seconds_per_day
       return sample_time.years * seconds_per_year + \
              sample_time.days * seconds_per_day + \
              sample_time.hours * seconds_per_hours + \
              sample_time.minutes * seconds_per_minutes + \
              sample_time.seconds

"""INPUTS"""

ell = 1
emm = 1
sample_time = SampleTime(0, 0, 1, 0, 0)

"""IMPORTS"""

import logging as log
log.basicConfig(level=log.INFO)
import toast
import toast.todmap
import toast.pipeline_tools
from toast.mpi import MPI
from toast_litebird.hardware import Hardware
from toast_litebird.focalplane import sim_telescope_detectors
from litebird_sim import Imo
import litebird_sim as lbs
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
log.info("python import succesfull")

"""ENVIRONMENT"""

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
    outdir = "./output/"
    sky_file = "./data/maptotlm_freq0.fits"
    beam_file = "./output/blm.fits"
# Focalplane = dictionary : [detector (key)] -> [quaternion tuple (value)]
focalplane = {list(hw.data['detectors'].items())[0][0]:list(hw.data['detectors'].items())[0][1]}
# Get detector list
detectors = sorted(focalplane.keys())
# Get quaterion tuple list
detquats = {}
for d in detectors:
    detquats[d] = focalplane[d]["quat"]

# One second of convolution
nsample = int(get_total_sample_time_seconds(sample_time) * args.sample_rate)
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

log.info("environment loaded")

"""COMPUTE"""

#itération sur les 15 fréquences:
for freq in range(1):

       log.info("computing for freq : %d" %freq)

       """PERTURBATION"""

       m = hp.read_map("./data/sim_sources_map_freq%s.fits"%freq)
       log.info("perturbation map loaded")

       hp.mollview(m, title = "Perturbation")
       plt.show()

       """SYNTHETIC BEAM"""

       data = toast.Data(comm)
       data.obs.append(obs) # 'append' is necessray (don't equal)

       nside = 64 
       lmax = nside # Pour retrouver : carte = lmax = nside requis
       nside_high = nside
       npix_high = 12 * nside_high ** 2
       lmax_high = nside_high * 2

       num_pix = hp.nside2npix(nside)
       ipix = np.arange(num_pix)

       raw_map_data = np.zeros(num_pix)
       array_alm_coeff = hp.map2alm(raw_map_data,lmax=lmax)

       # Change all values from the array of alm indexed by 'ell' and 'emm' previously found
       valchange1 = hp.Alm.getidx(lmax=lmax, l=ell, m=emm)
       array_alm_coeff[valchange1] = 1e-3

       #r = hp.rotator.Rotator(rot =(0,90,0))
       map_beam = hp.alm2map(array_alm_coeff,nside = nside,lmax=lmax)
       #map_beam = r.rotate_map_pixel(mapbeam)

       bl, blm = hp.anafast(np.array([map_beam,map_beam,map_beam]), lmax=lmax_high, iter=0, alm=True)
       hp.write_alm(args.beam_file, blm, overwrite=True)

       log.info("synthetic beam created")

       hp.mollview(map_beam,norm='hist', title = "Synthetic beam")
       plt.show()

       """CONVOLUTION"""

       toast.todmap.OpPointingHpix(nside=args.nside, nest=True, mode="IQU").exec(data)

       # Get HitMap (where the satellite sees during sample)
       num_pix = 12 * args.nside ** 2
       hitmap = np.zeros(num_pix)
       tod = data.obs[0]["tod"]
       for det in tod.local_dets:
           pixels = tod.cache.reference("pixels_{}".format(det))
           hitmap[pixels] = 1
       hitmap[hitmap == 0] = hp.UNSEEN

       # Show HitMap
       hp.mollview(hitmap, nest=True, cbar=False, title = "HitMap")
       hp.graticule(22.5, verbose=False)
       plt.show()

       name = "signal"
       toast.tod.OpCacheClear(name).exec(data) # clear 'name' entry to not overwrite
       conviqt = toast.todmap.OpSimConviqt(
           comm.comm_rank,
           args.sky_file,
           args.beam_file,
           lmax=2*nside_high,  # will use maximum from file
           beammmax=2*nside_high,  # will use maximum from file
           pol=True,
           fwhm=0, # width of a symmetric gaussian beam in arcmin
           order=13,
           calibrate=True, # normalize beam ?
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

       log.info("convolution performed")

       """WRITE DATA"""

       # From Tutorial (toast/tutorial/04_Simulated_Instrument_Signal/) :
       # "Destripe the signal and make a map. We use the nascent TOAST mapmaker because it can be run in serial mode without MPI. 
       # The TOAST mapmaker is still significantly slower so production runs should used libMadam.""
       mapmaker = toast.todmap.OpMapMaker(
           nside=args.nside,
           nnz=3,
           outdir=args.outdir,
           outprefix="toast_test_freq%s_ell%s_emm%s_"%(freq,ell,emm),
           baseline_length=10,
           iter_max=100,
           use_noise_prior=False,
       )
       mapmaker.exec(data)

       log.info("maps written")

       """PLOT FIGURES"""

       plt.figure(num = "Results")
       hits = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_hits.fits"%(freq,ell,emm))
       #hits[hits == 0] = hp.pixelfunc.UNSEEN
       hp.mollview(hits, sub=[2, 2, 1], title="hits")

       binmap = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_binned.fits"%(freq,ell,emm))
       #binmap[binmap == 0] = hp.pixelfunc.UNSEEN
       hp.mollview(binmap, sub=[2, 2, 2], title="binned map", cmap="coolwarm")

       destriped = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_destriped.fits"%(freq,ell,emm))
       #destriped[destriped == 0] = hp.pixelfunc.UNSEEN
       hp.mollview(destriped, sub=[2, 2, 3], title="destriped map", cmap="coolwarm")
       inmap = hp.ud_grade(hp.read_map("./data/sim_sources_map_freq%s.fits"%freq), args.nside)

       #inmap[hitmap == hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
       hp.mollview(inmap, sub=[2, 2, 4], title="input map", cmap="coolwarm")
       plt.savefig(args.outdir+"conviqt_freq%s_ell%s_emm%s.pdf"%(freq,ell,emm))
       plt.show()

log.info("compute finished")