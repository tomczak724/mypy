
import mypy
import pylab
from astropy.io import fits
from threedhst import eazyPy
from matplotlib.backends.backend_pdf import PdfPages


version = 'sg0023+0423_v0.1.5'
cat = mypy.readzout('/Users/atomczak/DATA/ORELSE/%s.zspec.cat' % version)
zout = mypy.readzout('/Users/atomczak/DATA/ORELSE/EAZY/%s.zout' % version)
fout = mypy.readcat('/Users/atomczak/DATA/ORELSE/%s.fout' % version)


imnames = ['B', 'V', 'R+', 'I+', 'r\'', 'i\'', 'z\'', 'J', 'K', '[3.6]', '[4.5]']
i_red, i_grn, i_blu = 3, 2, 1

imdats = [fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_B.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_V.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_Rsup.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_Isup.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_r.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_i.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_z.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_K.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/conv_noBackground_K.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_I1.fits'),
          fits.getdata('/Users/atomczak/DATA/ORELSE/wreg_Cl_0023+0423_I1.fits')]






pdf = PdfPages('galaxy-bios.pdf')




fig = pylab.figure(figsize=(10.89, 8.91))
sp_sed = pylab.subplot2grid((13, 8), (7, 0), rowspan=7, colspan=4)
sp_pz = pylab.subplot2grid((13, 8), (7, 5), rowspan=7, colspan=3)

sp_rgb = pylab.subplot2grid((13, 8), (0, 5), rowspan=6, colspan=3, xticks=[], yticks=[], aspect=1)

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

dr = 17
cutouts, bgs, noises = [], [], []
for i_im in range(len(imnames)):

	sp = sp_cutouts[i_im]

	subimdat = imdats[i_im][5400:5600, 5500:5700]
	bg = pylab.median(subimdat)
	noise = mypy.nmad(subimdat)
	bgs.append(bg)
	noises.append(noise)

	vmin = bg - 3 * noise
	vmax = bg + 25 * noise

	cutouts.append(sp.imshow(pylab.ones((2*dr+1, 2*dr+1)), interpolation='nearest', cmap=pylab.cm.gray_r, vmin=vmin, vmax=vmax, zorder=1))
	t = sp.text(0.03, 0.07, imnames[i_im], transform=sp.transAxes, color='r', fontweight='bold', zorder=2)


dummi = pylab.array([pylab.ones((2*dr+1, 2*dr+1)), pylab.ones((2*dr+1, 2*dr+1)), pylab.ones((2*dr+1, 2*dr+1))])
rgb = sp_rgb.imshow(dummi.transpose(), interpolation='nearest', cmap=pylab.cm.gray_r, zorder=1)




fig.subplots_adjust(hspace=0, wspace=0)




for i in [99]:


	###  thumbnail subplots
	xi, yi = int(cat.x[i] + 0.5), int(cat.y[i] + 0.5)
	for i_im in range(len(imnames)):
		cutouts[i_im].set_data(imdats[i_im][yi-dr:yi+dr+1, xi-dr:xi+dr+1])



	###  RGB subplot
	cutout_r = (imdats[i_red][yi-dr:yi+dr+1, xi-dr:xi+dr+1] - bgs[i_red]) / noises[i_red]
	cutout_g = (imdats[i_grn][yi-dr:yi+dr+1, xi-dr:xi+dr+1] - bgs[i_grn]) / noises[i_grn]
	cutout_b = (imdats[i_blu][yi-dr:yi+dr+1, xi-dr:xi+dr+1] - bgs[i_blu]) / noises[i_blu]
	stack = pylab.array([cutout_r, cutout_g, cutout_b])

	lo, hi = stack.min(), stack.max()
	stack -= lo
	stack /= (hi - lo)

	rgb.set_data(stack.transpose())




	###  sed subplot
	sed = eazyPy.getEazySED(i, MAIN_OUTPUT_FILE=version, OUTPUT_DIRECTORY='./EAZY')
	lambdaz, temp_sed, lci, temp_obs, fobs, efobs = sed

	sp_sed.loglog(lambdaz, temp_sed, color='gray', lw=1.5, zorder=1)
	sp_sed.errorbar(lci, temp_obs, ls='', marker='s', ms=12, mew=2, mec='gray', mfc='gray', zorder=2)
	sp_sed.errorbar(lci, fobs, yerr=efobs, ls='', marker='o', ms=12, mew=1.5, mfc='none', mec='r', ecolor='r', zorder=3)

	indspos = pylab.find(fobs > 0)
	fmin = min(fobs[indspos])
	fmax = max(fobs[indspos])

	sp_sed.set_xlim(2*10**3, 10**5)
	sp_sed.set_ylim(fmin/5, fmax*5)

	sp_sed.set_xlabel('wavelength')
	sp_sed.set_ylabel('f')




	###  P(z) subplot
	pz = eazyPy.getEazyPz(i, MAIN_OUTPUT_FILE=version, OUTPUT_DIRECTORY='./EAZY')

	sp_pz.fill_between(pz[0], 0, pz[1], color='#00b300')




















