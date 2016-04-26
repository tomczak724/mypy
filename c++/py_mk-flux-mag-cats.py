
import mypy
import pylab
from collections import namedtuple



field_name = 'sg0023+0423'
version = 'v0.2.1'

print '\nreading in catalogs...'
sexcats = [mypy.read_sexcat('../images/sextractor/catalog_B.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_V.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_Rsup.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_Isup.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_r.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_i.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_z.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_J.cat'),
           mypy.read_sexcat('../images/sextractor/catalog_K.cat')]

detect_i = 2
detection_sexcat = mypy.read_sexcat('../images/sextractor/catalog_detect_RIsup.cat')

iraccat = mypy.readcat('../images/tphot_test_RIsup-detect/combine.cat')
#irac1_cat = mypy.read_sexcat('../images/tphot_test_riz-detect/LRI_I1_SingleFit/lores_I1_tphot.cat_pass1_best')
#irac2_cat = mypy.read_sexcat('../images/tphot_test_riz-detect/LRI_I2_SingleFit/lores_I2_tphot.cat_pass1_best')

weight_ch1 = pylab.loadtxt('../images/sextractor/weights_I1.txt')
weight_ch2 = pylab.loadtxt('../images/sextractor/weights_I2.txt')

irac1_zp = 10**6 * (2.91e-6)**2 * 10**((25 - 8.9) / 2.5)
irac2_zp = 10**6 * (2.91e-6)**2 * 10**((25 - 8.9) / 2.5)

weightcats = [pylab.loadtxt('../images/sextractor/weights_B.txt'),
              pylab.loadtxt('../images/sextractor/weights_V.txt'),
              pylab.loadtxt('../images/sextractor/weights_Rsup.txt'),
              pylab.loadtxt('../images/sextractor/weights_Isup.txt'),
              pylab.loadtxt('../images/sextractor/weights_r.txt'),
              pylab.loadtxt('../images/sextractor/weights_i.txt'),
              pylab.loadtxt('../images/sextractor/weights_z.txt'),
              pylab.loadtxt('../images/sextractor/weights_J.txt'),
              pylab.loadtxt('../images/sextractor/weights_K.txt')]
print 'done!\n'


#  note: had to convert JK bands from Vega to AB due to 2MASS calibration
zeropoints = pylab.array([25., \
                          25., \
                          27.44, \
                          27.29, \
                          26.50, \
                          26.20, \
                          25.10, \
                          28.17 + 0.922, \
                          27.55 + 1.900])

#  zeropoint offsets derived from iterative template fitting
#zp_offsets = pylab.array([0.86691, 0.84541, 0.79648, 0.81173, 0.82525, 0.88421, 1.00])
zp_offsets = pylab.ones(len(zeropoints))
zeropoints -= 2.5 * pylab.log10(zp_offsets)



###  Background/Noise properties:
###    Fluxes were measured in ~10**4 randomly placed apertures
###    in each image with objects masked. The resulting histograms
###    were well-characterized by Gaussians and so the best-fit
###    mu and sigma are used as estimates of the background offset
###    and noise for the aperture size of the catalog.

### offsets and errors for r=6.5 pix apers
#bg_offsets = pylab.array([-0.0287, -0.0382, -0.5688, -0.6523, -0.3298, -0.3745, -0.3139, -15.41, -27.30])
bg_offsets = pylab.array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
flux_errors = pylab.array([0.0414, 0.0625, 1.0644, 1.5525, 0.5309, 0.6829, 0.9143, 108.63, 178.05])

### offsets and errors for r=3.25 pix apers
#bg_offsets = pylab.array([-0.0065581, -0.0158432, -0.3983771, -0.1517271, -0.0547522, -0.0598510, -0.0504124, -3.3135316, -6.3589188])
#bg_offsets = pylab.array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
#flux_errors = pylab.array([0.0158432, 0.0254563, 0.3983771, 0.6211592, 0.1768099, 0.2369858, 0.3129545, 47.4030118, 76.2378162])


#irac_offsets = pylab.array([-0.0207869, -0.0207869])
irac_offsets = pylab.array([0., 0.])
irac_fluxerr = pylab.array([0.0579433, 0.0745273])



MW_extinction_mag = pylab.array([0.078, 0.056, 0.049, 0.036, 0.049, 0.036, 0.027, 0.015, 0.006])
MW_extinction_scale = 10**(MW_extinction_mag / 2.5)







###  These are attempts at color-terms for LFC data as discussed with Roy
def get_rmag_true(rflux_instr, iflux_instr):
	a = 1.
	b = 0.05033836
	c = 26.52083003
	if rflux_instr < 0 or iflux_instr < 0:
		return -99.
	else:
		rmag_instr = -2.5 * pylab.log10(rflux_instr)
		imag_instr = -2.5 * pylab.log10(iflux_instr)
		return a * rmag_instr + b * (rmag_instr - imag_instr) + c

def get_imag_true(rflux_instr, iflux_instr):
	a = 1.
	b = 0.021948947
	c = 26.2291181
	if rflux_instr < 0 or iflux_instr < 0:
		return -99.
	else:
		rmag_instr = -2.5 * pylab.log10(rflux_instr)
		imag_instr = -2.5 * pylab.log10(iflux_instr)
		return a * imag_instr + b * (rmag_instr - imag_instr) + c

def get_zmag_true(iflux_instr, zflux_instr):
	a = 1.
	b = -0.03321189
	c = 25.02113495
	if iflux_instr < 0 or zflux_instr < 0:
		return -99.
	else:
		imag_instr = -2.5 * pylab.log10(iflux_instr)
		zmag_instr = -2.5 * pylab.log10(zflux_instr)
		return a * zmag_instr + b * (imag_instr - zmag_instr) + c


rmag_true_auto = pylab.array([get_rmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[4].flux_auto, sexcats[5].flux_auto)])
imag_true_auto = pylab.array([get_imag_true(f_r, f_i) for f_r, f_i in zip(sexcats[4].flux_auto, sexcats[5].flux_auto)])
zmag_true_auto = pylab.array([get_zmag_true(f_i, f_z) for f_i, f_z in zip(sexcats[5].flux_auto, sexcats[6].flux_auto)])

rmag_true_aper = pylab.array([get_rmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[4].flux_aper, sexcats[5].flux_aper)])
imag_true_aper = pylab.array([get_imag_true(f_r, f_i) for f_r, f_i in zip(sexcats[4].flux_aper, sexcats[5].flux_aper)])
zmag_true_aper = pylab.array([get_zmag_true(f_i, f_z) for f_i, f_z in zip(sexcats[5].flux_aper, sexcats[6].flux_aper)])



rflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(rmag_true_auto != -99)
rflux_true_auto[inds_pos] = 10**((25. - rmag_true_auto[inds_pos]) / 2.5)

iflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(imag_true_auto != -99)
iflux_true_auto[inds_pos] = 10**((25. - imag_true_auto[inds_pos]) / 2.5)

zflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(zmag_true_auto != -99)
zflux_true_auto[inds_pos] = 10**((25. - zmag_true_auto[inds_pos]) / 2.5)



rflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(rmag_true_aper != -99)
rflux_true_aper[inds_pos] = 10**((25. - rmag_true_aper[inds_pos]) / 2.5)

iflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(imag_true_aper != -99)
iflux_true_aper[inds_pos] = 10**((25. - imag_true_aper[inds_pos]) / 2.5)

zflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(zmag_true_aper != -99)
zflux_true_aper[inds_pos] = 10**((25. - zmag_true_aper[inds_pos]) / 2.5)


riz_mag_true_auto = [rmag_true_auto, imag_true_auto, zmag_true_auto]
riz_mag_true_aper = [rmag_true_aper, imag_true_aper, zmag_true_aper]
riz_flux_true_auto = [rflux_true_auto, iflux_true_auto, zflux_true_auto]
riz_flux_true_aper = [rflux_true_aper, iflux_true_aper, zflux_true_aper]




###  These are attempts at color-terms for R+I+ Subaru data as discussed with Roy
def get_Rplusmag_true(rflux_instr, iflux_instr):
	a = 1.
	b = 0.03185165
	c = 27.45145076
	if rflux_instr < 0 or iflux_instr < 0:
		return -99.
	else:
		rmag_instr = -2.5 * pylab.log10(rflux_instr)
		imag_instr = -2.5 * pylab.log10(iflux_instr)
		return a * rmag_instr + b * (rmag_instr - imag_instr) + c

def get_Iplusmag_true(rflux_instr, iflux_instr):
	a = 1.
	b = 0.00013198858
	c = 27.3378864
	if rflux_instr < 0 or iflux_instr < 0:
		return -99.
	else:
		rmag_instr = -2.5 * pylab.log10(rflux_instr)
		imag_instr = -2.5 * pylab.log10(iflux_instr)
		return a * imag_instr + b * (rmag_instr - imag_instr) + c


rplusmag_true_auto = pylab.array([get_Rplusmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[2].flux_auto, sexcats[3].flux_auto)])
iplusmag_true_auto = pylab.array([get_Iplusmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[2].flux_auto, sexcats[3].flux_auto)])

rplusmag_true_aper = pylab.array([get_Rplusmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[2].flux_aper, sexcats[3].flux_aper)])
iplusmag_true_aper = pylab.array([get_Iplusmag_true(f_r, f_i) for f_r, f_i in zip(sexcats[2].flux_aper, sexcats[3].flux_aper)])



rplusflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(rplusmag_true_auto != -99)
rplusflux_true_auto[inds_pos] = 10**((25. - rplusmag_true_auto[inds_pos]) / 2.5)

iplusflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(iplusmag_true_auto != -99)
iplusflux_true_auto[inds_pos] = 10**((25. - iplusmag_true_auto[inds_pos]) / 2.5)



rplusflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(rplusmag_true_aper != -99)
rplusflux_true_aper[inds_pos] = 10**((25. - rplusmag_true_aper[inds_pos]) / 2.5)

iplusflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(iplusmag_true_aper != -99)
iplusflux_true_aper[inds_pos] = 10**((25. - iplusmag_true_aper[inds_pos]) / 2.5)


riplus_mag_true_auto = [rplusmag_true_auto, iplusmag_true_auto]
riplus_mag_true_aper = [rplusmag_true_aper, iplusmag_true_aper]
riplus_flux_true_auto = [rplusflux_true_auto, iplusflux_true_auto]
riplus_flux_true_aper = [rplusflux_true_aper, iplusflux_true_aper]



###  These are attempts at color-terms for JK wfcam data as discussed with Roy
###  The parameterizations below are already converted to the AB system
def get_Jmag_true(jflux_instr, kflux_instr):
	a = 1.
	b = -0.04157139
	c = 29.07118567
	if jflux_instr < 0 or kflux_instr < 0:
		return -99.
	else:
		jmag_instr = -2.5 * pylab.log10(jflux_instr)
		kmag_instr = -2.5 * pylab.log10(kflux_instr)
		return a * jmag_instr + b * (jmag_instr - kmag_instr) + c

def get_Kmag_true(jflux_instr, kflux_instr):
	a = 1.
	b = -0.0257698323
	c = 29.4494709
	if jflux_instr < 0 or kflux_instr < 0:
		return -99.
	else:
		jmag_instr = -2.5 * pylab.log10(jflux_instr)
		kmag_instr = -2.5 * pylab.log10(kflux_instr)
		return a * kmag_instr + b * (jmag_instr - kmag_instr) + c


jmag_true_auto = pylab.array([get_Jmag_true(f_j, f_k) for f_j, f_k in zip(sexcats[7].flux_auto, sexcats[8].flux_auto)])
kmag_true_auto = pylab.array([get_Kmag_true(f_j, f_k) for f_j, f_k in zip(sexcats[7].flux_auto, sexcats[8].flux_auto)])

jmag_true_aper = pylab.array([get_Jmag_true(f_j, f_k) for f_j, f_k in zip(sexcats[7].flux_aper, sexcats[8].flux_aper)])
kmag_true_aper = pylab.array([get_Kmag_true(f_j, f_k) for f_j, f_k in zip(sexcats[7].flux_aper, sexcats[8].flux_aper)])



jflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(jmag_true_auto != -99)
jflux_true_auto[inds_pos] = 10**((25. - jmag_true_auto[inds_pos]) / 2.5)

kflux_true_auto = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(kmag_true_auto != -99)
kflux_true_auto[inds_pos] = 10**((25. - kmag_true_auto[inds_pos]) / 2.5)



jflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(jmag_true_aper != -99)
jflux_true_aper[inds_pos] = 10**((25. - jmag_true_aper[inds_pos]) / 2.5)

kflux_true_aper = -99. * pylab.ones(len(detection_sexcat.number))
inds_pos = pylab.find(kmag_true_aper != -99)
kflux_true_aper[inds_pos] = 10**((25. - kmag_true_aper[inds_pos]) / 2.5)


jk_mag_true_auto = [jmag_true_auto, kmag_true_auto]
jk_mag_true_aper = [jmag_true_aper, kmag_true_aper]
jk_flux_true_auto = [jflux_true_auto, kflux_true_auto]
jk_flux_true_aper = [jflux_true_aper, kflux_true_aper]








translate_name = field_name + '_' + version + '.translate'
print 'generating', translate_name, '\n'

outer_translate = open(translate_name, 'w')
outer_translate.write('fluxaper_B   F8\nerraper_B   E8\n')
outer_translate.write('fluxaper_V   F9\nerraper_V   E9\n')
outer_translate.write('fluxaper_Rplus  F14\nerraper_Rplus  E14\n')
outer_translate.write('fluxauto_Rplus  TOT14\n')
outer_translate.write('fluxaper_Iplus  F15\nerraper_Iplus  E15\n')
outer_translate.write('fluxaper_r   F5\nerraper_r   E5\n')
outer_translate.write('fluxaper_i   F6\nerraper_i   E6\n')
outer_translate.write('fluxaper_z   F7\nerraper_z   E7\n')
outer_translate.write('fluxaper_J   F17\nerraper_J   E17\n')
outer_translate.write('fluxaper_K   F19\nerraper_K   E19\n')
outer_translate.write('fluxcolor_ch1  F20\nerrcolor_ch1  E20\n')
outer_translate.write('fluxcolor_ch2  F21\nerrcolor_ch2  E21\n')
outer_translate.close()





###########################################################
###  Adding z_specs from B. Lemaux's sepctroscopic catalog
###########################################################

specfile = '../external_catalogs/FINAL.SG0023.deimos.lris.feb2012.nodups.cat'

colnames = ['id', 'slitmask', 'slitnum', 'ra', 'dec', 'mag1', 'mag2', 'mag3', 'z_spec', 'z_spec_err', 'Q']
speccat_obj = namedtuple('speccat_obj', colnames)
speccat = speccat_obj(*pylab.loadtxt(specfile, dtype=str, usecols=range(len(colnames)), unpack=1))

id_specs = (pylab.zeros(len(sexcats[0].number)) - 1).tolist()
z_specs = pylab.zeros(len(sexcats[0].number)) - 1
z_spec_errs = pylab.zeros(len(sexcats[0].number)) - 1
Q_specs = pylab.zeros(len(sexcats[0].number)) - 99

print 'searching for spectroscopic matches...'
inds_match_tomczak = []
inds_match_lemaux = []
nmatch = []



for i in range(len(speccat.id)):
	
	###  Skipping serendipitous objects with no previously identified photometric counterpart
	if i > 0 and speccat.ra[i] == speccat.ra[i-1] and speccat.dec[i] == speccat.dec[i-1]:
		continue

	ra, dec = float(speccat.ra[i]), float(speccat.dec[i])
	seps = mypy.radec_sep(ra, dec, sexcats[0].x_world, sexcats[0].y_world)

	inds = pylab.find(seps < 0.75)
	nmatch.append(len(inds))

	if len(inds) == 1:
		id_specs[inds[0]] = speccat.id[i]
		z_specs[inds[0]] = float(speccat.z_spec[i])
		z_spec_errs[inds[0]] = float(speccat.z_spec_err[i])
		Q_specs[inds[0]] = int(speccat.Q[i])
		inds_match_tomczak.append(inds[0])
		inds_match_lemaux.append(i)

	if len(inds) == 0:
		print speccat.ra[i], speccat.dec[i], speccat.Q[i], speccat.slitnum[i], speccat.id[i]


inds_match_tomczak = pylab.array(inds_match_tomczak)
inds_match_lemaux = pylab.array(inds_match_lemaux)
nmatch = pylab.array(nmatch)
print 'done!\n'









outer_match = open('%s_%s.crossmatch' % (field_name, version), 'w')
outer_match.write('# id_spec slitmask slitnum ra_spec dec_spec z_spec z_spec_err Q id_phot ra_phot dec_phot nmatch blend\n')


for i in range(len(speccat.id)):

	blendflag = 0
	if i > 0 and speccat.ra[i] == speccat.ra[i-1] and speccat.dec[i] == speccat.dec[i-1]:
		blendflag = 1


	ra, dec = float(speccat.ra[i]), float(speccat.dec[i])
	seps = mypy.radec_sep(ra, dec, sexcats[0].x_world, sexcats[0].y_world)

	inds = pylab.find(seps < 0.75)

	if len(inds) == 0:
		id_phot = '%7i' % -1
		ra_phot = '-1.000000'
		dec_phot = '-1.000000'
		n = 0

	if len(inds) == 1:
		id_phot = '%7i' % sexcats[0].number[inds[0]]
		ra_phot = '%.7f' % sexcats[0].x_world[inds[0]]
		dec_phot = '%.7f' % sexcats[0].y_world[inds[0]]
		n = 1

	if len(inds) > 1:
		n = len(inds)
		blendflag = 1
		id_phot = '%i' % sexcats[0].number[inds[0]]
		ra_phot = '%.7f' % sexcats[0].x_world[inds[0]]
		dec_phot = '%.7f' % sexcats[0].y_world[inds[0]]
		for ind in inds[1:]:
			id_phot += ',%i' % sexcats[0].number[ind]
			ra_phot += ',%.7f' % sexcats[0].x_world[ind]
			dec_phot += ',%.7f' % sexcats[0].y_world[ind]

	outer_match.write('%-13s' % speccat.id[i])
	outer_match.write(' %-8s' % speccat.slitmask[i])
	outer_match.write(' %-7s' % speccat.slitnum[i])
	outer_match.write(' %11s' % speccat.ra[i])
	outer_match.write(' %11s' % speccat.dec[i])
	outer_match.write(' %10s' % speccat.z_spec[i])
	outer_match.write(' %11s' % speccat.z_spec_err[i])
	outer_match.write(' %3s' % speccat.Q[i])

	outer_match.write(' %s' % id_phot)
	outer_match.write('  %s' % ra_phot)
	outer_match.write('  %s' % dec_phot)
	outer_match.write(' %2i' % n)
	outer_match.write(' %2i' % blendflag)
	outer_match.write('\n')

outer_match.close()








#       Catalog Format
#  id
#  x
#  y
#  ra
#  dec
#  fluxaper_B, erraper_B, fluxaper_V, erraper_V ... fluxaper_K, erraper_K                       # aperture fluxes
#  ftot_B, etot_B, ftot_V, etot_V ... ftot_K, etot_K     # total fluxes
#  w_B, w_V ... w_K                                      # weights
#  wmin                                                  # minimum weight


#  *.cat
outer_flux = open(field_name + '_' + version + '.cat', 'w')
outer_flux.write('# id z_spec x y ra dec ')

outer_flux.write('fluxaper_B erraper_B ')
outer_flux.write('fluxaper_V erraper_V ')
outer_flux.write('fluxaper_Rplus erraper_Rplus ')
outer_flux.write('fluxaper_Iplus erraper_Iplus ')
outer_flux.write('fluxaper_r erraper_r ')
outer_flux.write('fluxaper_i erraper_i ')
outer_flux.write('fluxaper_z erraper_z ')
outer_flux.write('fluxaper_J erraper_J ')
outer_flux.write('fluxaper_K erraper_K ')
outer_flux.write('fluxcolor_ch1 errcolor_ch1 ')
outer_flux.write('fluxcolor_ch2 errcolor_ch2 ')

outer_flux.write('fluxauto_B errauto_B ')
outer_flux.write('fluxauto_V errauto_V ')
outer_flux.write('fluxauto_Rplus errauto_Rplus ')
outer_flux.write('fluxauto_Iplus errauto_Iplus ')
outer_flux.write('fluxauto_r errauto_r ')
outer_flux.write('fluxauto_i errauto_i ')
outer_flux.write('fluxauto_z errauto_z ')
outer_flux.write('fluxauto_J errauto_J ')
outer_flux.write('fluxauto_K errauto_K ')

outer_flux.write('weight_B weight_V weight_Rplus weight_Iplus weight_r weight_i weight_z weight_J weight_K weight_ch1 weight_ch2 wmin')
outer_flux.write('\n')



#  *.mag
outer_mag = open(field_name + '_' + version + '.mag', 'w')
outer_mag.write('# id z_spec ra dec ')

outer_mag.write('magaper_B erraper_B ')
outer_mag.write('magaper_V erraper_V ')
outer_mag.write('magaper_Rplus erraper_Rplus ')
outer_mag.write('magaper_Iplus erraper_Iplus ')
outer_mag.write('magaper_r erraper_r ')
outer_mag.write('magaper_i erraper_i ')
outer_mag.write('magaper_z erraper_z ')
outer_mag.write('magaper_J erraper_J ')
outer_mag.write('magaper_K erraper_K ')
outer_mag.write('magcolor_ch1 errcolor_ch1 ')
outer_mag.write('magcolor_ch2 errcolor_ch2 ')

outer_mag.write('apercorr ')

outer_mag.write('weight_B weight_V weight_Rplus weight_Iplus weight_r weight_i weight_z weight_J weight_K weight_ch1 weight_ch2 wmin')
outer_mag.write('\n')



#  *.zspec.cat
outer_flux_zspec = open(field_name + '_' + version + '.zspec.cat', 'w')
outer_flux_zspec.write('# id z_spec x y ra dec ')

outer_flux_zspec.write('fluxaper_B erraper_B ')
outer_flux_zspec.write('fluxaper_V erraper_V ')
outer_flux_zspec.write('fluxaper_Rplus erraper_Rplus ')
outer_flux_zspec.write('fluxaper_Iplus erraper_Iplus ')
outer_flux_zspec.write('fluxaper_r erraper_r ')
outer_flux_zspec.write('fluxaper_i erraper_i ')
outer_flux_zspec.write('fluxaper_z erraper_z ')
outer_flux_zspec.write('fluxaper_J erraper_J ')
outer_flux_zspec.write('fluxaper_K erraper_K ')
outer_flux_zspec.write('fluxcolor_ch1 errcolor_ch1 ')
outer_flux_zspec.write('fluxcolor_ch2 errcolor_ch2 ')

outer_flux_zspec.write('fluxauto_B errauto_B ')
outer_flux_zspec.write('fluxauto_V errauto_V ')
outer_flux_zspec.write('fluxauto_Rplus errauto_Rplus ')
outer_flux_zspec.write('fluxauto_Iplus errauto_Iplus ')
outer_flux_zspec.write('fluxauto_r errauto_r ')
outer_flux_zspec.write('fluxauto_i errauto_i ')
outer_flux_zspec.write('fluxauto_z errauto_z ')
outer_flux_zspec.write('fluxauto_J errauto_J ')
outer_flux_zspec.write('fluxauto_K errauto_K ')

outer_flux_zspec.write('weight_B weight_V weight_Rplus weight_Iplus weight_r weight_i weight_z weight_J weight_K weight_ch1 weight_ch2 wmin')
outer_flux_zspec.write('\n')



#  *.match
outer_match = open(field_name + '_' + version + '.match', 'w')
outer_match.write('# id_phot id_spec ra dec z_spec z_spec_err Q_spec')
outer_match.write('\n')


#  looping over number of objects
#print '\nDISCLAIMER: scaling IRAC total fluxes to aperture by r-band total to aperture ratio\n'

print 'generating flux catalog...'
for object_i in range(len(sexcats[0].number)):

	####################################
	###  0th pass for the match catalog
	####################################
	outer_match.write('% 7i' % sexcats[0].number[object_i])
	outer_match.write(' %12s' % id_specs[object_i])

	outer_match.write(' % 3.7f' % sexcats[0].x_world[object_i])
	outer_match.write(' % 3.7f' % sexcats[0].y_world[object_i])

	outer_match.write(' %.7f' % z_specs[object_i])
	outer_match.write(' %.7f' % z_spec_errs[object_i])
	outer_match.write(' %3i' % Q_specs[object_i])
	outer_match.write('\n')




	##################################
	###  First pass for flux catalog
	##################################
	outer_flux.write('% 7i' % sexcats[0].number[object_i])

	if Q_specs[object_i] in [3, 4]:
		outer_flux.write(' %.5f' % z_specs[object_i])
	else:
		outer_flux.write(' %.5f' % -1)

	outer_flux.write(' % 6.2f' % sexcats[0].x_image[object_i])
	outer_flux.write(' % 6.2f' % sexcats[0].y_image[object_i])

	outer_flux.write(' % 3.7f' % sexcats[0].x_world[object_i])
	outer_flux.write(' % 3.7f' % sexcats[0].y_world[object_i])



	#  looping over filters
	#  first, aperture fluxes and errors
	for filter_i in range(len(sexcats)):

		###  write -99s for objects with weight < 0.05 of median
		if weightcats[filter_i][object_i] < 0.05:
			outer_flux.write(' -99.0000')
			outer_flux.write(' -99.0000')
		else:
			sexcat = sexcats[filter_i]
			scale_ab25 = 10**((25 - zeropoints[filter_i]) / 2.5)

			### grab LFC color-termed fluxes
			if filter_i in [4, 5, 6]:
				flux_i = riz_flux_true_aper[filter_i-4][object_i]
			### ... or Subaru R+I+ color-termed fluxes
			elif filter_i in [2, 3]:
				flux_i = riplus_flux_true_aper[filter_i-2][object_i]
			### ... or WFCAM JK color-termed fluxes
			elif filter_i in [7, 8]:
				flux_i = jk_flux_true_aper[filter_i-7][object_i]
			else:
				flux_i = sexcat.flux_aper[object_i]
				flux_i -= bg_offsets[filter_i]
				flux_i *= scale_ab25
			flux_i *= MW_extinction_scale[filter_i]

			eflux_i = flux_errors[filter_i]
			eflux_i /= weightcats[filter_i][object_i]**0.5
			eflux_i *= scale_ab25
			eflux_i *= MW_extinction_scale[filter_i]

			outer_flux.write(' %.7f' % flux_i)
			outer_flux.write(' %.7f' % eflux_i)


	# adding IRAC fluxes
	if detection_sexcat.flux_auto[object_i] > 0 and detection_sexcat.flux_aper[object_i] > 0:
		scale_color = detection_sexcat.flux_auto[object_i] / detection_sexcat.flux_aper[object_i]
	else:
		scale_color = 1.

	flux_i = iraccat.ch1[object_i]
	flux_i -= irac_offsets[0]
	flux_i *= irac1_zp
	flux_i /= scale_color    # it may be the case that this introduces undue scatter

	eflux_i = irac_fluxerr[0] * 1.
	eflux_i /= weight_ch1[object_i]**0.5
	eflux_i *= irac1_zp
	eflux_i /= scale_color    # it may be the case that this introduces undue scatter

	if weight_ch1[object_i] < 0.05:
		outer_flux.write(' -99.0000')
		outer_flux.write(' -99.0000')
	else:
		outer_flux.write(' %.7f' % flux_i)
		outer_flux.write(' %.7f' % eflux_i)


	flux_i = iraccat.ch2[object_i]
	flux_i -= irac_offsets[1]
	flux_i *= irac2_zp
	flux_i /= scale_color    # it may be the case that this introduces undue scatter

	eflux_i = irac_fluxerr[1] * 1.
	eflux_i /= weight_ch2[object_i]**0.5
	eflux_i *= irac2_zp
	eflux_i /= scale_color    # it may be the case that this introduces undue scatter

	if weight_ch2[object_i] < 0.05:
		outer_flux.write(' -99.0000')
		outer_flux.write(' -99.0000')
	else:
		outer_flux.write(' %.7f' % flux_i)
		outer_flux.write(' %.7f' % eflux_i)




	#  second, total fluxes and errors
	for filter_i in range(len(sexcats)):

		###  write -99s for objects with weight < 0.05 of median
		if weightcats[filter_i][object_i] < 0.05:
			outer_flux.write(' -99.0000')
			outer_flux.write(' -99.0000')
		else:
			sexcat = sexcats[filter_i]
			scale_ab25 = 10**((25 - zeropoints[filter_i]) / 2.5)

			### grab LFC color-termed fluxes
			if filter_i in [4, 5, 6]:
				flux_i = riz_flux_true_auto[filter_i-4][object_i]
			### ... or Subaru R+I+ color-termed fluxes
			elif filter_i in [2, 3]:
				flux_i = riplus_flux_true_auto[filter_i-2][object_i]
			### ... or WFCAM JK color-termed fluxes
			elif filter_i in [7, 8]:
				flux_i = jk_flux_true_auto[filter_i-7][object_i]
			else:
				flux_i = sexcat.flux_auto[object_i]
				flux_i -= bg_offsets[filter_i]
				flux_i *= scale_ab25
			flux_i *= MW_extinction_scale[filter_i]

			eflux_i = flux_errors[filter_i]
			eflux_i /= weightcats[filter_i][object_i]**0.5
			eflux_i *= scale_ab25
			eflux_i *= MW_extinction_scale[filter_i]

			outer_flux.write(' %.7f' % flux_i)
			outer_flux.write(' %.7f' % eflux_i)


	#  last, weights
	weights = []
	for filter_i in range(len(sexcats)):
		weightcat = weightcats[filter_i]
		outer_flux.write(' %.2f' % weightcat[object_i])
		weights.append(weightcat[object_i])

	outer_flux.write(' %.2f' % weight_ch1[object_i])
	outer_flux.write(' %.2f' % weight_ch2[object_i])
	weights.append(weight_ch1[object_i])
	weights.append(weight_ch2[object_i])

	outer_flux.write(' %.2f' % min(weights))
	outer_flux.write('\n')







	##################################
	#  Second pass for magnitude catalog
	##################################
	outer_mag.write('% 7i' % sexcats[0].number[object_i])

	if Q_specs[object_i] in [3, 4]:
		outer_mag.write(' %.5f' % z_specs[object_i])
	else:
		outer_mag.write(' %.5f' % -1)

	outer_mag.write(' % 3.7f' % sexcats[0].x_world[object_i])
	outer_mag.write(' % 3.7f' % sexcats[0].y_world[object_i])


	#  looping over filters
	#  first, aperture fluxes and errors
	for filter_i in range(len(sexcats)):
		sexcat = sexcats[filter_i]
		scale_ab25 = 10**((25 - zeropoints[filter_i]) / 2.5)

		### grab LFC color-termed fluxes
		if filter_i in [4, 5, 6]:
			flux_i = riz_flux_true_aper[filter_i-4][object_i]
		### ... or Subaru R+I+ color-termed fluxes
		elif filter_i in [2, 3]:
			flux_i = riplus_flux_true_aper[filter_i-2][object_i]
		### ... or WFCAM JK color-termed fluxes
		elif filter_i in [7, 8]:
			flux_i = jk_flux_true_aper[filter_i-7][object_i]
		else:
			flux_i = sexcat.flux_auto[object_i]
			flux_i -= bg_offsets[filter_i]
			flux_i *= scale_ab25
		flux_i *= MW_extinction_scale[filter_i]

		eflux_i = flux_errors[filter_i] / weightcats[filter_i][object_i]**0.5 * scale_ab25 * MW_extinction_scale[filter_i]
		weight_i = weightcats[filter_i][object_i]

		if flux_i < 0:
			mag_i_s = ' -99.00'
			emag_i_s = ' -99.0'
		else:
			mag_i = -2.5 * pylab.log10(flux_i) + 25
			mag_i_s = ' %.3f' % mag_i
			if weight_i > 0:
				emag_i = eflux_i * 2.5 / flux_i / pylab.log(10)
				emag_i_s = ' %.3f' % emag_i
			else:
				emag_i_s = ' -99.0'

		outer_mag.write(mag_i_s)
		outer_mag.write(emag_i_s)

	#  calculating aperture correction
	if detection_sexcat.flux_auto[object_i] > 0 and detection_sexcat.flux_aper[object_i] > 0:
		scale_color = detection_sexcat.flux_auto[object_i] / detection_sexcat.flux_aper[object_i]
		apercorr = -2.5 * pylab.log10(detection_sexcat.flux_auto[object_i] / detection_sexcat.flux_aper[object_i])
		apercorr_s = ' %.3f' % apercorr
	else:
		scale_color = 1.
		apercorr_s = ' 0.0000'


	# adding IRAC fluxes
	flux_i = iraccat.ch1[object_i]
	flux_i -= irac_offsets[0]
	flux_i *= irac1_zp
	flux_i /= scale_color    # it may be the case that this introduces undue scatter

	eflux_i = irac_fluxerr[0] * 1.
	eflux_i /= weight_ch1[object_i]**0.5
	eflux_i *= irac1_zp
	eflux_i /= scale_color    # it may be the case that this introduces undue scatter

	if flux_i < 0:
		mag_i_s = ' -99.00'
		emag_i_s = ' -99.0'
	else:
		mag_i = -2.5 * pylab.log10(flux_i) + 25
		mag_i_s = ' %.3f' % mag_i
		if weight_i > 0:
			emag_i = eflux_i * 2.5 / flux_i / pylab.log(10)
			emag_i_s = ' %.3f' % emag_i
		else:
			emag_i_s = ' -99.0'

	outer_mag.write(mag_i_s)
	outer_mag.write(emag_i_s)



	flux_i = iraccat.ch2[object_i]
	flux_i -= irac_offsets[1]
	flux_i *= irac2_zp
	flux_i /= scale_color    # it may be the case that this introduces undue scatter

	eflux_i = irac_fluxerr[1] * 1.
	eflux_i /= weight_ch2[object_i]**0.5
	eflux_i *= irac2_zp
	eflux_i /= scale_color    # it may be the case that this introduces undue scatter

	if flux_i < 0:
		mag_i_s = ' -99.00'
		emag_i_s = ' -99.0'
	else:
		mag_i = -2.5 * pylab.log10(flux_i) + 25
		mag_i_s = ' %.3f' % mag_i
		if weight_i > 0:
			emag_i = eflux_i * 2.5 / flux_i / pylab.log(10)
			emag_i_s = ' %.3f' % emag_i
		else:
			emag_i_s = ' -99.0'

	outer_mag.write(mag_i_s)
	outer_mag.write(emag_i_s)

	outer_mag.write(apercorr_s)



	#  last, weights
	weights = []
	for filter_i in range(len(sexcats)):
		weightcat = weightcats[filter_i]
		outer_mag.write(' %.2f' % weightcat[object_i])
		weights.append(weightcat[object_i])

	outer_mag.write(' %.2f' % weight_ch1[object_i])
	outer_mag.write(' %.2f' % weight_ch2[object_i])
	weights.append(weight_ch1[object_i])
	weights.append(weight_ch2[object_i])

	outer_mag.write(' %.2f' % min(weights))
	outer_mag.write('\n')







	##################################
	#  Lastly for the z_spec only catalog
	##################################
	if Q_specs[object_i] in [3, 4]:
		outer_flux_zspec.write('% 7i' % sexcats[0].number[object_i])

		outer_flux_zspec.write(' %.5f' % z_specs[object_i])

		outer_flux_zspec.write(' % 6.2f' % sexcats[0].x_image[object_i])
		outer_flux_zspec.write(' % 6.2f' % sexcats[0].y_image[object_i])
		outer_flux_zspec.write(' % 3.7f' % sexcats[0].x_world[object_i])
		outer_flux_zspec.write(' % 3.7f' % sexcats[0].y_world[object_i])

		#  looping over filters
		#  first, aperture fluxes and errors
		for filter_i in range(len(sexcats)):
		
		###  write -99s for objects with weight < 0.05 of median
			if weightcats[filter_i][object_i] < 0.05:
				outer_flux_zspec.write(' -99.0000')
				outer_flux_zspec.write(' -99.0000')
			else:
				sexcat = sexcats[filter_i]
				scale_ab25 = 10**((25 - zeropoints[filter_i]) / 2.5)

				### grab LFC color-termed fluxes
				if filter_i in [4, 5, 6]:
					flux_i = riz_flux_true_aper[filter_i-4][object_i]
				### ... or Subaru R+I+ color-termed fluxes
				elif filter_i in [2, 3]:
					flux_i = riplus_flux_true_aper[filter_i-2][object_i]
				### ... or WFCAM JK color-termed fluxes
				elif filter_i in [7, 8]:
					flux_i = jk_flux_true_aper[filter_i-7][object_i]
				else:
					flux_i = sexcat.flux_aper[object_i]
					flux_i -= bg_offsets[filter_i]
					flux_i *= scale_ab25
				flux_i *= MW_extinction_scale[filter_i]

				eflux_i = flux_errors[filter_i]
				eflux_i /= weightcats[filter_i][object_i]**0.5
				eflux_i *= scale_ab25
				eflux_i *= MW_extinction_scale[filter_i]

				outer_flux_zspec.write(' %.7f' % flux_i)
				outer_flux_zspec.write(' %.7f' % eflux_i)


		# adding IRAC fluxes
		if detection_sexcat.flux_auto[object_i] > 0 and detection_sexcat.flux_aper[object_i] > 0:
			scale_color = detection_sexcat.flux_auto[object_i] / detection_sexcat.flux_aper[object_i]
		else:
			scale_color = 1.

		flux_i = iraccat.ch1[object_i]
		flux_i -= irac_offsets[0]
		flux_i *= irac1_zp
		flux_i /= scale_color    # it may be the case that this introduces undue scatter

		eflux_i = irac_fluxerr[0] * 1.
		eflux_i /= weight_ch1[object_i]**0.5
		eflux_i *= irac1_zp
		eflux_i /= scale_color    # it may be the case that this introduces undue scatter

		if weight_ch1[object_i] < 0.05:
			outer_flux.write(' -99.0000')
			outer_flux.write(' -99.0000')
		else:
			outer_flux_zspec.write(' %.7f' % flux_i)
			outer_flux_zspec.write(' %.7f' % eflux_i)


		flux_i = iraccat.ch2[object_i]
		flux_i -= irac_offsets[1]
		flux_i *= irac2_zp
		flux_i /= scale_color    # it may be the case that this introduces undue scatter

		eflux_i = irac_fluxerr[1] * 1.
		eflux_i /= weight_ch2[object_i]**0.5
		eflux_i *= irac2_zp
		eflux_i /= scale_color    # it may be the case that this introduces undue scatter

		if weight_ch2[object_i] < 0.05:
			outer_flux.write(' -99.0000')
			outer_flux.write(' -99.0000')
		else:
			outer_flux_zspec.write(' %.7f' % flux_i)
			outer_flux_zspec.write(' %.7f' % eflux_i)


		#  second, total fluxes and errors
		for filter_i in range(len(sexcats)):

		###  write -99s for objects with weight < 0.05 of median
			if weightcats[filter_i][object_i] < 0.05:
				outer_flux_zspec.write(' -99.0000')
				outer_flux_zspec.write(' -99.0000')
			else:
				sexcat = sexcats[filter_i]
				scale_ab25 = 10**((25 - zeropoints[filter_i]) / 2.5)

				### grab LFC color-termed fluxes
				if filter_i in [4, 5, 6]:
					flux_i = riz_flux_true_auto[filter_i-4][object_i]
				### ... or Subaru R+I+ color-termed fluxes
				elif filter_i in [2, 3]:
					flux_i = riplus_flux_true_auto[filter_i-2][object_i]
				### ... or WFCAM JK color-termed fluxes
				elif filter_i in [7, 8]:
					flux_i = jk_flux_true_auto[filter_i-7][object_i]
				else:
					flux_i = sexcat.flux_auto[object_i]
					flux_i -= bg_offsets[filter_i]
					flux_i *= scale_ab25
				flux_i *= MW_extinction_scale[filter_i]

				eflux_i = flux_errors[filter_i]
				eflux_i /= weightcats[filter_i][object_i]**0.5
				eflux_i *= scale_ab25
				eflux_i *= MW_extinction_scale[filter_i]

				outer_flux_zspec.write(' %.7f' % flux_i)
				outer_flux_zspec.write(' %.7f' % eflux_i)


		#  last, weights
		weights = []
		for filter_i in range(len(sexcats)):
			weightcat = weightcats[filter_i]
			outer_flux_zspec.write(' %.2f' % weightcat[object_i])
			weights.append(weightcat[object_i])

		outer_flux_zspec.write(' %.2f' % weight_ch1[object_i])
		outer_flux_zspec.write(' %.2f' % weight_ch2[object_i])
		weights.append(weight_ch1[object_i])
		weights.append(weight_ch2[object_i])

		outer_flux_zspec.write(' %.2f' % min(weights))
		outer_flux_zspec.write('\n')



outer_match.close()
outer_flux.close()
outer_mag.close()
outer_flux_zspec.close()
print 'done!\n'
print 'wrote to: ', field_name + '_' + version + '.cat'
print 'wrote to: ', field_name + '_' + version + '.mag'
print 'wrote to: ', field_name + '_' + version + '.zspec.cat'
print 'wrote to: ', field_name + '_' + version + '.match\n'




















