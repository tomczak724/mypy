
import os
import sys
import time
import read
import dusts
import ezgal
import numpy
import extras
import pickle
import filters
import gen_chi2_cube
from astropy import cosmology

t_start = time.time()



if len(sys.argv) > 1: paramfile = sys.argv[1]
else: paramfile = 'fast.param'


###  reading parameter file
params = read.readparams(paramfile)
print '\nRead parameter file: %s' % paramfile

###  setting cosmology
cosmo = cosmology.FlatLambdaCDM(H0=float(params.H0), Om0=float(params.OMEGA_M))
z_cosmo = numpy.arange(0., 10., 0.01)
time_cosmo = cosmo.age(z_cosmo).value


###  reading filter files
res = read.readres(params.FILTERS_RES)
translate = numpy.loadtxt('%s.translate' % params.CATALOG, dtype=str, comments='#')
print 'Read filter files: %s, %s.translate' % (params.FILTERS_RES, params.CATALOG)


###  reading flux catalog
flux_catalog = numpy.genfromtxt('%s.cat' % params.CATALOG, names=True)
print 'Read flux catalog: %s.cat' % params.CATALOG


###  reading zout catalog
zout_catalog = numpy.genfromtxt('%s.zout' % params.CATALOG, names=True)
print 'Read EAZY catalog: %s.zout' % params.CATALOG





###  generating model cube/gridspace
grid_z = numpy.arange(float(params.Z_MIN), float(params.Z_MAX), float(params.Z_STEP))
grid_lage = numpy.arange(float(params.LOG_AGE_MIN), float(params.LOG_AGE_MAX), float(params.LOG_AGE_STEP))
grid_ltau = numpy.arange(float(params.LOG_TAU_MIN), float(params.LOG_TAU_MAX), float(params.LOG_TAU_STEP))
grid_Av = numpy.arange(float(params.A_V_MIN), float(params.A_V_MAX), float(params.A_V_STEP))

n_z = len(grid_z)
n_lage = len(grid_lage)
n_ltau = len(grid_ltau)
n_Av = len(grid_Av)



#############################
###  reading in model files
#############################
filter_lams = []
filter_files = []
model_files = []
path2models = '%s/ised_%s.%s' % (params.LIBRARY_DIR, params.SFH, params.RESOLUTION)
for i_ltau in xrange(len(grid_ltau)):

	ltau = grid_ltau[i_ltau]
	modelname = '%s/%s_%s_%s_z%s_ltau%.1f.ised' % (path2models, params.LIBRARY, params.RESOLUTION, params.IMF, params.METAL[2:], ltau)
	modelobj = ezgal.model(modelname)
	modelobj.set_cosmology(Om=float(params.OMEGA_M), Ol=float(params.OMEGA_L), h=float(params.H0) / 100., w=-1)

	#ages_

	###  writing temporary filter files
	for i_filt in xrange(len(translate)):
		if translate[i_filt][1][0] == 'F':
			i_res_filt = int(translate[i_filt][1][1:]) - 1
			res_filt = res[i_res_filt]
			numpy.savetxt(translate[i_filt][1], res_filt.data)

			###  adding filter to model
			modelobj.add_filter(str(translate[i_filt][1]))
			os.remove(translate[i_filt][1])

			if i_ltau == 0:
				res_filt.name = translate[i_filt][0]
				filter_files.append(res_filt)
				filter_lams.append(filters.get_filter_central_wavelength(res_filt.lams, res_filt.trans))

		###  adding name of error column
		if translate[i_filt][1][0] == 'E':
			if i_ltau == 0:
				res_filt.err_name = translate[i_filt][0]



	modelobj.zs = grid_z
	model_files.append(modelobj)

filter_lams = numpy.array(filter_lams)
n_filters = len(filter_files)





#####################################
###  generating template error cube
#####################################
template_err_file = numpy.loadtxt(params.TEMP_ERR_FILE)
template_err_cube = numpy.zeros((n_z, n_filters))

for i_z in range(n_z):
	t_zlam = template_err_file[:,0] * (1 + grid_z[i_z])
	t_err = template_err_file[:,1]
	template_err_file_zi = numpy.array(zip(t_zlam, t_err))
	for i_filt in range(n_filters):
		template_err_cube[i_z][i_filt] += filters.synphot(template_err_file_zi, filter_files[i_filt].data)






############################
###  generating model cube
############################
class model_cube_class:
	def __init__(self, grid_z, grid_lage, grid_ltau, grid_Av):
		self.grid_z = grid_z
		self.grid_lage = grid_lage
		self.grid_ltau = grid_ltau
		self.grid_Av = grid_Av
		self.flux_cube = numpy.zeros((n_z, n_ltau, n_lage, n_Av, n_filters))
		self.mass_cube = numpy.zeros((n_ltau, n_lage, n_Av))
		self.sfh_cube = numpy.zeros((n_ltau, n_lage, n_Av))

try:
	model_cube = pickle.load(open('model_cube.pickle', 'rb'))
except:
	model_cube = model_cube_class(grid_z, grid_lage, grid_ltau, grid_Av)
	for i_z in xrange(n_z):
		print i_z, n_z
		for i_ltau in xrange(n_ltau):
			model = model_files[i_ltau]

			for i_lage in xrange(n_lage):
				sed_wavelengths = model.ls
				sed_flux = model.get_sed(10**grid_lage[i_lage] / 10**9., age_units='gyrs', units='Fv')

				for i_Av in xrange(n_Av):
					Alam2Av = dusts.calzetti(sed_wavelengths)
					corr = 10**(-0.4 * grid_Av[i_Av] * Alam2Av)
					sed_flux_corr = sed_flux * corr
					sed = numpy.array(zip(sed_wavelengths[::-1] * (1 + grid_z[i_z]), sed_flux_corr[::-1]))

					for i_filt in xrange(n_filters):
						model_cube.flux_cube[i_z][i_ltau][i_lage][i_Av][i_filt] += filters.synphot(sed, filter_files[i_filt].data) * extras.lum2flux(grid_z[i_z], cosmo, ZP_AB=float(params.AB_ZEROPOINT))

					model_cube.mass_cube[i_ltau][i_lage][i_Av] = numpy.interp(10**grid_lage[i_lage], model.ages, model.masses)
					model_cube.sfh_cube[i_ltau][i_lage][i_Av] = numpy.interp(10**grid_lage[i_lage], model.ages, model.sfh)

	pickle.dump(model_cube, open('model_cube.pickle', 'wb'))









#############################################
###  testing writing code to calculate chi2
#############################################

i_gals = numpy.array([ 1029,  1127,  1244,  1318,  1384,  1510,  1549,  1771,  1836,
        1935,  2022,  2052,  2057,  2620,  2743,  2784,  2867,  2868,
        3147,  3178,  3225,  3253,  3291,  3347,  3430,  3449,  3464,
        3498,  3530,  3534,  3600,  3631,  3651,  3655,  3660,  3665,
        3668,  3690,  3708,  3731,  3808,  3812,  3818,  3850,  3851,
        3857,  3936,  3946,  3947,  3998,  4015,  4026,  4029,  4052,
        4057,  4079,  4080,  4094,  4136,  4169,  4214,  4307,  4311,
        4317,  4321,  4350,  4354,  4390,  4393,  4438,  4497,  4572,
        4603,  4711,  4762,  4783,  4831,  4859,  4936,  4958,  5099,
        5114,  5120,  5156,  5240,  5297,  5305,  5313,  5373,  5384,
        5572,  5581,  5591,  5597,  5607,  5662,  5713,  5725,  5744,
        5762,  5803,  5840,  5891,  5899,  5908,  5931,  5938,  5960,
        6018,  6052,  6059,  6079,  6096,  6097,  6122,  6149,  6223,
        6286,  6288,  6294,  6340,  6348,  6361,  6393,  6406,  6418,
        6481,  6486,  6488,  6496,  6506,  6513,  6520,  6546,  6604,
        6624,  6729,  6740,  6774,  6794,  6810,  6811,  6849,  6853,
        6863,  6875,  6932,  6978,  7004,  7046,  7093,  7095,  7123,
        7152,  7169,  7172,  7203,  7204,  7277,  7313,  7323,  7338,
        7410,  7422,  7447,  7470,  7496,  7548,  7563,  7566,  7573,
        7574,  7614,  7625,  7635,  7666,  7675,  7676,  7750,  7764,
        7838,  7845,  7887,  7917,  7935,  8011,  8042,  8108,  8113,
        8169,  8189,  8215,  8222,  8230,  8299,  8343,  8355,  8433,
        8486,  8585,  8614,  8709,  8734,  8780,  8799,  8821,  8866,
        8875,  8886,  8906,  8910,  8935,  8941,  8979,  8995,  9025,
        9083,  9100,  9131,  9186,  9192,  9196,  9214,  9230,  9241,
        9273,  9328,  9340,  9366,  9404,  9422,  9546,  9574,  9578,
        9602,  9672,  9683,  9684,  9692,  9698,  9716,  9738,  9780,
        9795,  9809,  9815,  9821,  9831,  9895,  9922,  9943,  9953,
        9963,  9969, 10031, 10068, 10178, 10237, 10243, 10246, 10261,
       10267, 10289, 10369, 10411, 10422, 10463, 10491, 10494, 10505,
       10507, 10520, 10524, 10538, 10630, 10631, 10641, 10649, 10691,
       10721, 10735, 10752, 10787, 10794, 10814, 10868, 10890, 10936,
       10945, 10987, 11008, 11032, 11059, 11087, 11099, 11103, 11114,
       11115, 11150, 11165, 11210, 11284, 11341, 11355, 11358, 11401,
       11402, 11420, 11423, 11445, 11456, 11492, 11520, 11534, 11538,
       11558, 11561, 11588, 11601, 11620, 11625, 11650, 11673, 11678,
       11742, 11757, 11776, 11792, 11814, 11852, 11880, 11885, 11899,
       11924, 11963, 12010, 12028, 12029, 12037, 12045, 12077, 12078,
       12089, 12094, 12106, 12125, 12144, 12161, 12169, 12195, 12267,
       12270, 12275, 12314, 12329, 12345, 12348, 12351, 12357, 12367,
       12376, 12391, 12393, 12469, 12473, 12502, 12512, 12513, 12514,
       12557, 12573, 12588, 12607, 12613, 12651, 12688, 12718, 12744,
       12751, 12759, 12783, 12789, 12833, 12866, 12877, 12919, 12921,
       12946, 13024, 13056, 13064, 13116, 13172, 13180, 13198, 13234,
       13262, 13263, 13270, 13300, 13302, 13335, 13406, 13411, 13440,
       13492, 13506, 13508, 13521, 13527, 13587, 13589, 13616, 13640,
       13652, 13690, 13701, 13749, 13783, 13805, 13823, 13824, 13978,
       14014, 14043, 14059, 14211, 14306, 14309, 14323, 14450, 14481,
       14564, 14570, 14622, 14645, 14661, 14669, 14697, 14698, 14768,
       14798, 14898, 14903, 14937, 15010, 15042, 15088, 15164, 15187,
       15244, 15384, 15430, 15496, 15537, 15541, 15581, 15626, 15678,
       15948, 15959, 16070, 16230, 16237, 16295, 16441, 16443, 16769,
       16850, 16868, 17140, 17167, 17226, 17299, 17368, 17546, 17624,
       17940, 17994, 18002, 18099, 18134, 18480, 18595, 18908, 19085,
       19310, 19441, 19801, 20042, 20121, 20124, 20137, 20271])

bestfits = []
chi2red_cube = numpy.zeros((n_ltau, n_lage, n_Av))
scale_cube = numpy.zeros((n_ltau, n_lage, n_Av))


t0_p = time.time()

#for asdf, id_gal in enumerate(flux_catalog['id']):
#	i_gal = id_gal - 1
#	mypy.progress_bar(asdf, len(flux_catalog['id']))

for asdf, i_gal in enumerate(i_gals):

#for i_gal in [1127]:

	###  finding nearest redshift grid-point to galaxy's redshift
	if zout_catalog['z_spec'][i_gal] > 0:
		dz = grid_z - zout_catalog['z_spec'][i_gal]
	else:
		dz = grid_z - zout_catalog[params.NAME_ZPHOT][i_gal]

	i_z = numpy.argmin(abs(dz))



	fluxes = numpy.zeros(n_filters)
	efluxes = numpy.zeros(n_filters)
	for i_filt in xrange(n_filters):
		fluxes[i_filt] = flux_catalog[filter_files[i_filt].name][i_gal]
		efluxes[i_filt] = flux_catalog[filter_files[i_filt].err_name][i_gal]

	etot = numpy.sqrt(efluxes**2 + (template_err_cube[i_z]*fluxes)**2)
	weights = 1. / etot**2



	###  calculating scale factors between data and models
	for i_ltau in xrange(n_ltau):
		for i_lage in xrange(n_lage):
			for i_Av in xrange(n_Av):
				m = model_cube.flux_cube[i_z][i_ltau][i_lage][i_Av]

				wmm = numpy.sum(weights * m * m)
				wfm = numpy.sum(weights * fluxes * m)
				scale = wfm / wmm
				scale_cube[i_ltau][i_lage][i_Av] = scale
				chi2red_cube[i_ltau][i_lage][i_Av] = numpy.sum((fluxes - m*scale)**2 / etot**2 / n_filters)

	ibest_ltau, ibest_lage, ibest_Av = numpy.where(chi2red_cube == chi2red_cube.min())
	ibest_ltau, ibest_lage, ibest_Av = ibest_ltau[0], ibest_lage[0], ibest_Av[0]
	bestfits.append([ibest_ltau, ibest_lage, ibest_Av])
	#bestfits.append([grid_ltau[ibest_ltau], grid_lage[ibest_lage], grid_Av[ibest_Av]])

bestfits = numpy.array(bestfits)

tf_p = time.time()











t0_c = time.time()

for asdf, id_gal in enumerate(flux_catalog['id']):
	i_gal = id_gal - 1


	###  finding nearest redshift grid-point to galaxy's redshift
	if zout_catalog['z_spec'][i_gal] > 0:
		dz = grid_z - zout_catalog['z_spec'][i_gal]
	else:
		dz = grid_z - zout_catalog[params.NAME_ZPHOT][i_gal]

	i_z = numpy.argmin(abs(dz))



	fluxes = numpy.zeros(n_filters)
	efluxes = numpy.zeros(n_filters)
	for i_filt in xrange(n_filters):
		fluxes[i_filt] = flux_catalog[filter_files[i_filt].name][i_gal]
		efluxes[i_filt] = flux_catalog[filter_files[i_filt].err_name][i_gal]


	asdf = gen_chi2_cube.gen_chi2_cube(fluxes, efluxes, template_err_cube[i_z], model_cube.flux_cube[i_z])


tf_c = time.time()













fig = pylab.figure()
sp = fig.add_subplot(111)
sp.axis([2*10**3, 10**5, 10**-2, 2*10**2])
sp.set_xscale('log')
sp.set_yscale('log')

pylab.errorbar(filter_lams, fluxes, efluxes, ls='', marker='o')
modpoints = pylab.plot(filter_lams, m*scale, 'ro', mec='r', mew=3, ms=10, mfc='none')[0]



for iz in xrange(n_Av):
	for iy in xrange(n_ltau):
		for ix in xrange(n_lage):

			m = model_cube.flux_cube[i_z][iy][ix][iz]

			wmm = numpy.sum(weights * m * m)
			wfm = numpy.sum(weights * fluxes * m)
			scale = wfm / wmm

			modpoints.set_ydata(m*scale)
			pylab.draw()

			time.sleep(0.2)























