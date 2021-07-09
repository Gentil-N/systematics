
import pymaster as nmt 
import pysm
from pysm.nominal import models
import time

nside = 256
Npix = hp.nside2npix(nside)
Nf=9
#lmax = nside*3-1
lmax=850
scale = 5
Nlbin = 10
dusttype = 0

# Initialize binning scheme with Nlbin ells per bandpower
b = nmt.bins.NmtBin(nside=nside,lmax=lmax,nlb=Nlbin)
leff = b.get_effective_ells()

#carte
###################################################################################################################
#à utiliser une fois pour créer cartes de poussières differentes et les sauvegarder

sky_config = {
    'dust' : models("d%s"%dusttype, nside),
    'synchrotron' : models("s%s"%dusttype, nside),

#    'cmb' : models("c1",nside),
}

sky = pysm.Sky(sky_config)

litebird = np.load(Pr+"Librairies/instrus/litebird.npy",allow_pickle = True).item()
litebird.update({"nside" : nside})
litebird.update({"use_smoothing" : False})

instrument = pysm.Instrument(litebird)

# objet de taille [2,15,3,3000] : [propre/bruit, fréquences, IQU,pixels]

freq_maps = np.array(instrument.observe(sky, write_outputs=False))

FM=np.array(freq_maps)
map1=FM[0,:,:,:]
map2 = FM[1,:,:,:]
maptot = map1

cl, alm = hp.anafast(maptot[0], lmax=2*nside, iter=0, alm=True)
for f in range(len(maptot[:,0,0])):
    hp.write_alm(Pr+"systés/beam_synth+conviqt/conviqt_perturb/"+"maptotlm_freq%s.fits"%f, alm, overwrite=True)
    hp.fitsfunc.write_map(Pr+"systés/beam_synth+conviqt/conviqt_perturb/"+"sim_sources_map_freq%s.fits"%f,maptot[f], overwrite=True)
