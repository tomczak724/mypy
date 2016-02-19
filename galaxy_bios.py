
import mypy
import pylab
from astropy.io import fits
from threedhst import eazyPy
import matplotlib.patheffects as PathEffects
from matplotlib.backends.backend_pdf import PdfPages


def restore_logaxes_labels(subplot, xaxis=True, yaxis=True):
	'''
	This script reformats the default labeling scheme
	for logarithmic axis labels to the "regular" format.
	e.g. axis labels of "10**-1" and "10**3" will be
	changed to "0.1" and "1000" respectively.
	'''
	xticks = subplot.get_xticks() 
	yticks = subplot.get_yticks() 

	if xaxis:
		xticks_new = []
		for xi in xticks:
			base = pylab.log10(xi)
			if base >= 0:
				xi_new = '%i' % xi
			else:
				formatter = '%.' + str(int(abs(base))) + 'f'
				xi_new = formatter % xi
			xticks_new.append(xi_new)
		subplot.xaxis.set_ticklabels(xticks_new)

	if yaxis:
		yticks_new = []
		for yi in yticks:
			base = pylab.log10(yi)
			if base >= 0:
				yi_new = '%i' % yi
			else:
				formatter = '%.' + str(int(abs(base))) + 'f'
				yi_new = formatter % yi
			yticks_new.append(yi_new)
		subplot.yaxis.set_ticklabels(yticks_new)

	return subplot




def add_inset(subplot, rect=[0.4, 0.4, 0.2, 0.2]):
	'''
	This script creates an Axes instance within the
	coordinate frame of the provided subplot.
	'''

	###  coordinates of the subplot in the figure window's coordinate frame
	box = subplot.get_position()
	xlo, xhi = box.x0, box.x1
	ylo, yhi = box.y0, box.y1
	dx, dy = (xhi - xlo), (yhi - ylo)

	###  width/height of requested axes in the figure window's coordinate frame
	sub_dx = dx * rect[2]
	sub_dy = dy * rect[3]

	###  position of requested axes in the figure window's coordinate frame
	sub_xlo = xlo + dx * rect[0]
	sub_ylo = ylo + dy * rect[1]

	inset = subplot.figure.add_axes([sub_xlo, sub_ylo, sub_dx, sub_dy])
	return inset










version = 'sg0023+0423_v0.1.5'
cat = mypy.readzout('/Users/atomczak/DATA/ORELSE/%s.zspec.cat' % version)
zout = mypy.readzout('/Users/atomczak/DATA/ORELSE/EAZY/%s.zout' % version)
fout = mypy.readcat('/Users/atomczak/DATA/ORELSE/%s.fout' % version)


imnames = ['B', 'V', 'R$_+$', 'I$_+$', 'r\'', 'i\'', 'z\'', 'J', 'K', '[3.6]', '[4.5]']
i_red, i_grn, i_blu = 3, 2, 1

imdats = [fits.getdata('/Users/atomczak/DATA/ORELSE/Cl_0023+0423_tomczak_B_zp25.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/Cl_0023+0423_tomczak_V_zp25.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_CL0023_tomczak_R.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_CL0023_tomczak_I.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_cl0023_r_shift_norm_unwarp.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_cl0023_i_shift_norm_unwarp.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_cl0023_z_shift_norm_unwarp.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_J.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_K.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_I1.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_I2.fits')]







fig = pylab.figure(figsize=(10.89, 8.91))
sp_info = pylab.subplot2grid((13, 8), (7, 5), rowspan=7, colspan=5, xticks=[], yticks=[])

sp_sed = pylab.subplot2grid((13, 8), (7, 0), rowspan=7, colspan=5)

#sp_pz = pylab.subplot2grid((13, 8), (7, 5), rowspan=7, colspan=3)
sp_pz = add_inset(sp_sed, [0.7, 0.74, 0.29, 0.25])
sp_pz.set_yticks([])
sp_pz.xaxis.set_tick_params(labelsize=14)


sp_rgb = pylab.subplot2grid((13, 8), (0, 5), rowspan=6, colspan=3, xticks=[], yticks=[], aspect=1)
sp_rgb.text(0.03, 0.03, '%s %s %s' % (imnames[i_blu], imnames[i_grn], imnames[i_red]), \
	        size=20, transform=sp_rgb.transAxes, color='w', fontweight='bold', \
	        path_effects=[PathEffects.withStroke(linewidth=3, foreground='k')])

sp_cutouts = [pylab.subplot2grid((13, 8), (0, 0), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (0, 1), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (0, 2), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (0, 3), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (2, 0), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (2, 1), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (2, 2), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (2, 3), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (4, 0), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (4, 1), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1),
              pylab.subplot2grid((13, 8), (4, 2), rowspan=2, colspan=1, xticks=[], yticks=[], aspect=1)]


bg_V, noise_V = (-0.0007834, 0.003343)
bg_R, noise_R = (-0.0099281, 0.045694)
bg_I, noise_I = (-0.0103544, 0.067142)

pxscale = 0.2
dr_asec = 7. / 2
dr_px = int(dr_asec / pxscale)
cutouts, bgs, noises = [], [], []
for i_im in range(len(imnames)):

	sp = sp_cutouts[i_im]

	subimdat = imdats[i_im][5400:5600, 5500:5700]
	bg = pylab.median(subimdat)
	noise = mypy.nmad(subimdat)
	bgs.append(bg)
	noises.append(noise)

	vmin = bg - 3 * noise
	vmax = bg + 20 * noise

	cutouts.append(sp.imshow(pylab.ones((2*dr_px+1, 2*dr_px+1)), interpolation='nearest', cmap=pylab.cm.gray_r, vmin=vmin, vmax=vmax, zorder=1))
	ax = sp.axis()
	sp.axis([ax[0], ax[1], ax[3], ax[2]])

	t = sp.text(0.03, 0.07, imnames[i_im], transform=sp.transAxes, \
		        color='k', fontweight='bold', zorder=2, \
		        path_effects=[PathEffects.withStroke(linewidth=2, foreground='w')])


dr_px_rgb = 2 * dr_px
dummi = pylab.ones((2*dr_px_rgb+1, 2*dr_px_rgb+1, 3))
rgb = sp_rgb.imshow(dummi, interpolation='nearest', cmap=pylab.cm.gray_r, zorder=1)
ax = sp_rgb.axis()
sp_rgb.axis([ax[0], ax[1], ax[3], ax[2]])

sidelen_asec = dummi.shape[0] * pxscale
two_asec = 2. / sidelen_asec

sp_rgb.text(0.9-two_asec*2/3., 0.07, '2"', size=16, transform=sp_rgb.transAxes, color='w', fontweight='bold', \
	        path_effects=[PathEffects.withStroke(linewidth=2, foreground='k')])

sp_rgb.plot([0.9, 0.9-two_asec], [0.05, 0.05], color='k', lw=4, transform=sp_rgb.transAxes)
sp_rgb.plot([0.9, 0.9-two_asec], [0.05, 0.05], color='w', lw=2, transform=sp_rgb.transAxes)




fig.subplots_adjust(hspace=0, wspace=0)




for i in range(4):

	###  info subplot
	sp_info.text(0.05, 0.84, '%s\nid = %6i' % (version, cat.id[i]))

	use = 1
	if use: sp_info.text(0.05, 0.75, 'use = %i' % use)
	else: sp_info.text(0.05, 0.75, 'use = %i' % use, color='r', fontweight='bold')

	star = 0
	if not star: sp_info.text(0.05, 0.67, 'star = %i' % star)
	else: sp_info.text(0.05, 0.67, 'star = %i' % star, color='r', fontweight='bold')

	zp, dzlo, dzhi = zout.z_peak[i], zout.z_peak[i]-zout.l68[i], zout.u68[i]-zout.z_peak[i], 
	sp_info.text(0.05, 0.32, 'z$_{phot}$ = $%.2f^{\it{+}%.2f}_{\it{-}%.2f}$' % (zp, dzlo, dzhi) + \
	                        '\nz$_{spec}$ = %.4f       Q = %i' % (zout.z_spec[i], 3) + \
		                    '\nlog M$_*$ = %.2f' % fout.lmass[i])

	rflux = cat.fluxaper_Rplus[i]
	iflux = cat.fluxaper_Iplus[i]
	if rflux > 0: rmag = '%.1f' % (-2.5 * pl.log10(rflux) + 25)
	else: rmag = '-99'
	if iflux > 0: imag = '%.1f' % (-2.5 * pl.log10(iflux) + 25)
	else: imag = '-99'
	sp_info.text(0.05, 0.07, u'R$_+$ = %s\nI$_+$ = %s' % (rmag, imag))




	###  thumbnail subplots
	xi, yi = int(cat.x[i] + 0.5), int(cat.y[i] + 0.5)
	for i_im in range(len(imnames)):
		cutouts[i_im].set_data(imdats[i_im][yi-dr_px:yi+dr_px+1, xi-dr_px:xi+dr_px+1])



	###  RGB subplot
	cutout_r = (imdats[i_red][yi-dr_px_rgb:yi+dr_px_rgb+1, xi-dr_px_rgb:xi+dr_px_rgb+1] - bg_I) / noise_I
	cutout_g = (imdats[i_grn][yi-dr_px_rgb:yi+dr_px_rgb+1, xi-dr_px_rgb:xi+dr_px_rgb+1] - bg_R) / noise_R
	cutout_b = (imdats[i_blu][yi-dr_px_rgb:yi+dr_px_rgb+1, xi-dr_px_rgb:xi+dr_px_rgb+1] - bg_V) / noise_V
	stack = pylab.zeros((cutout_r.shape[1], cutout_r.shape[0], 3))
	stack[:,:,0] = cutout_r
	stack[:,:,1] = cutout_g
	stack[:,:,2] = cutout_b

	lo, hi = stack.min(), stack.max()
	stack -= lo
	stack /= (hi - lo)

	stack *= 2
	hinds = pylab.where(stack > 1)
	stack[hinds] = 1

	rgb.set_data(stack)




	###  sed subplot
	sed = eazyPy.getEazySED(i, MAIN_OUTPUT_FILE=version, OUTPUT_DIRECTORY='./EAZY')
	lambdaz, temp_sed, lci, temp_obs, fobs, efobs = sed
	lambdaz, lci = lambdaz / 10.**4, lci / 10.**4

	sp_sed.loglog(lambdaz, temp_sed, color='r', lw=1, zorder=1)
#	sp_sed.errorbar(lci, temp_obs, ls='', marker='s', ms=8, mew=2, mec='r', mfc='r', zorder=2)
	sp_sed.errorbar(lci, fobs, yerr=efobs, ls='', marker='o', ms=12, mew=2, mfc='none', mec='k', ecolor='k', zorder=3)

	indspos = pylab.find(fobs > 0)
	fmin = min(fobs[indspos])
	fmax = max(fobs[indspos])

	sp_sed.set_xlim(2*10**-1, 1.7*10**1)
	sp_sed.set_ylim(fmin/5, fmax*5)
	sp_sed = restore_logaxes_labels(sp_sed)

	sp_sed.set_xlabel('wavelength [$\mu m$]')
	sp_sed.set_ylabel('F lambda')




	###  P(z) subplot
	sp_pz.set_xlabel('z', size=13)
	sp_pz.set_ylabel('P(z)', rotation=90, size=13)
	
	pz = eazyPy.getEazyPz(i, MAIN_OUTPUT_FILE=version, OUTPUT_DIRECTORY='./EAZY')

	sp_pz.fill_between(pz[0], 0, pz[1], color='#00b300')
	sp_pz.set_xlim(0, 5)




	pylab.savefig('/Users/atomczak/mugshot%i_%s.pdf' % (cat.id[i], version))



	sp_info.clear()
	sp_sed.clear()
	sp_pz.clear()

	sp_info.set_xticks([])
	sp_info.set_yticks([])
	sp_pz.set_yticks([])









