"""INPUTS"""

ell = 1
emm = 1

"""IMPORTS"""

# Capture C++ output in the jupyter cells
#%reload_ext wurlitzer
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

mpiworld, procs, rank = toast.mpi.get_world()
comm = toast.mpi.Comm(mpiworld)
# Appel de l'imo pour avoir les infos sur les detecteurs et la scanning strategy: 
imo = Imo(flatfile_location="./data")
scanning_strategy = lbs.SpinningScanningStrategy.from_imo(imo, "5e3bf49f-53f3-4f80-a65c-6191a4820b86")
#appel d'un fichier qui contient les detecteurs: (je n'en ai pris qu'une, idéalement il faudrait toutes les prendre)
hw = Hardware("./data/hardware_config_files/20Hardware.toml.gz")
log.info("detector loaded")

class args:
    sample_rate = 19  # Hz
    hwp_rpm = None
    hwp_step_deg = None
    hwp_step_time_s = None
    spin_period_min = 60*scanning_strategy.spin_rate_hz # 10
    spin_angle_deg = 50
    prec_period_min = 60*scanning_strategy.precession_rate_hz # 50
    prec_angle_deg = 45 
    coord = "E"
    nside = 256
    nnz = 3
    outdir = "./output"
    sky_file = outdir + "maptotlm.fits"
    beam_file = outdir + "blm.fits"
#On crée une focaleplane avec les quaternions des detectors
focalplane = {list(hw.data['detectors'].items())[0][0]:list(hw.data['detectors'].items())[0][1]}
#focalplane = fake_focalplane(samplerate=args.sample_rate, fknee=0.1, alpha=2)
detectors = sorted(focalplane.keys())
detquats = {}
for d in detectors:
    detquats[d] = focalplane[d]["quat"]
#On convolue sur un an:
nsample = int(1*365.25*24*60*60*args.sample_rate)
log.info("focalplane created")

start_sample = 0
start_time = 0
iobs = 0

tod = toast.todmap.TODSatellite(
        comm.comm_group,
        detquats,
        nsample,
        coord=args.coord,
        firstsamp=start_sample,
        firsttime=start_time,
        rate=args.sample_rate,
        spinperiod=args.spin_period_min,
        spinangle=args.spin_angle_deg,
        precperiod=args.prec_period_min,
        precangle=args.prec_angle_deg,
        detranks=comm.group_size,
        hwprpm=args.hwp_rpm,
        hwpstep=args.hwp_step_deg,
        hwpsteptime=args.hwp_step_time_s,
    )

# Constantly slewing precession axis                                                                                                                                             
precquat = np.empty(4 * tod.local_samples[1], dtype=np.float64).reshape((-1, 4))
toast.todmap.slew_precession_axis(
        precquat,
        firstsamp=start_sample + tod.local_samples[0],
        samplerate=args.sample_rate,
        degday=360.0 / 365.25,)
tod.set_prec_axis(qprec=precquat)
#noise = toast.pipeline_tools.get_analytic_noise(args, comm, focalplane)
obs = {}
obs["name"] = "science_{:05d}".format(iobs)
obs["tod"] = tod
obs["intervals"] = None
obs["baselines"] = None
#obs["noise"] = noise
obs["id"] = iobs
# Conviqt requires at least minimal focal plane information to be present in the observation
obs["focalplane"] = toast.pipeline_tools.Focalplane(focalplane)
log.info("observation's data created")

log.info("environment loaded")

"""COMPUTE"""

#itération sur les 15 fréquences:
for freq in range(15):

       log.info("computing for freq : %d" %freq)

       """PERTURBATION"""

       #J'appelle la carte qu'on va convoluer avec de la poussière et du synchrotron dedans
       data = toast.Data(comm)
       data.obs.append(obs)
       nside_high = 256
       npix_high = 12 * nside_high ** 2
       lmax_high = nside_high * 2
       m = hp.read_map(args.outdir+"sim_sources_map_freq%s.fits"%freq)
       # hp.mollview(m)
       # plt.show()

       log.info("perturbation map loaded")

       """SYNTHETIC BEAM"""

       #Je crée les beams perturbés:

       nside = 256 
       lmax = nside #Pour retrouver carte = lmax = nside requis 

       npix = hp.nside2npix(nside)
       ipix = np.arange(npix)

       map0 = np.zeros(npix)
       almtest = hp.map2alm(map0,lmax=lmax)

       #L,M = hp.Alm.getlm(lmax=lmax)
       # ell = np.arange(lmax+1)
       # L = []
       # M = []
       # for i in range(len(ell)):
       #     M.append(list(np.arange(ell[i]+1)))
       #     L.append(ell[i]*np.ones(len(np.arange(ell[i]+1))))
       # L = np.array([item for sublist in L for item in sublist]).astype(int)
       # M = np.array([item for sublist in M for item in sublist])
       # lenlmax= len(almtest)
       #numb = hp.Alm.getidx(lmax=lmax, l=L, m=M)

       almtest2 = 0*almtest 
       valchange1 = hp.Alm.getidx(lmax=lmax, l=ell, m=emm)
       almtest2[valchange1] = 1e-3

       r = hp.rotator.Rotator(rot =(0,90,0))
       mapbeam = hp.alm2map(almtest2,nside = nside,lmax=lmax)
       mapbeam = r.rotate_map_pixel(mapbeam)

       #décomente ça si tu veux voir le beam
       # hp.mollview(mapbeam,norm='hist')
       # plt.show()

       bl, blm = hp.anafast(np.array([mapbeam,mapbeam,mapbeam]), lmax=lmax_high, iter=0, alm=True)
       hp.write_alm(args.beam_file, blm, overwrite=True)

       log.info("synthetic beam created")

       """CONVOLUTION"""

       #J'appel conviqt exactement comme dans toast/tutorial/04_Simulated_Instrument_Signal/")"

       toast.todmap.OpPointingHpix(nside=args.nside, nest=True, mode="IQU").exec(data)

       npix = 12 * args.nside ** 2
       hitmap = np.zeros(npix)
       tod = data.obs[0]["tod"]
       for det in tod.local_dets:
           pixels = tod.cache.reference("pixels_{}".format(det))
           hitmap[pixels] = 1
       hitmap[hitmap == 0] = hp.UNSEEN
       #hp.mollview(hitmap, nest=True, title="all hit pixels", cbar=False)
       #plt.show()
       hp.graticule(22.5, verbose=False)

       name = "signal"
       toast.tod.OpCacheClear(name).exec(data)
       conviqt = toast.todmap.OpSimConviqt(
           comm.comm_rank,
           args.sky_file,
           args.beam_file,
           lmax=2*nside_high,  # Will use maximum from file
           beammmax=2*nside_high,  # Will use maximum from file
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

       log.info("convolution performed")

       """WRITE DATA"""

       mapmaker = toast.todmap.OpMapMaker(
           nside=args.nside,
           nnz=3,
           #name=name,
           outdir=args.outdir,
           outprefix="toast_test_freq%s_ell%s_emm%s_"%(freq,ell,emm),
           baseline_length=10,
           # maskfile=self.maskfile_binary,
           # weightmapfile=self.maskfile_smooth,
           # subharmonic_order=None,
           iter_max=100,
           use_noise_prior=False,
           # precond_width=30,
       )
       mapmaker.exec(data)

       log.info("map written")

       plt.figure(figsize=[12, 8])

       hitmap = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_hits.fits"%(freq,ell,emm))
       hitmap[hitmap == 0] = hp.UNSEEN
       #hp.mollview(hitmap, sub=[2, 2, 1], title="hits")

       binmap = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_binned.fits"%(freq,ell,emm))
       binmap[binmap == 0] = hp.UNSEEN
       #hp.mollview(binmap, sub=[2, 2, 2], title="binned map", cmap="coolwarm")

       destriped = hp.read_map(args.outdir+"/toast_test_freq%s_ell%s_emm%s_destriped.fits"%(freq,ell,emm))
       destriped[destriped == 0] = hp.UNSEEN
       #hp.mollview(destriped, sub=[2, 2, 3], title="destriped map", cmap="coolwarm")
       inmap = hp.ud_grade(hp.read_map(args.outdir+"sim_sources_map.fits"), args.nside)

       inmap[hitmap == hp.UNSEEN] = hp.UNSEEN
       #hp.mollview(inmap, sub=[2, 2, 4], title="input map", cmap="coolwarm")
       plt.savefig(args.outdir+"conviqt_freq%s_ell%s_emm%s.pdf"%(freq,ell,emm))
       #plt.show()

log.info("compute finished")