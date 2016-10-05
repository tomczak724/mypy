
###  This script is designed to take an image and estimate
###  a limiting magnitude by inserting fake sources and 
###  attempting to recover them using the same SE parameters
###  as for the images. The fake sources are a 2D sersic
###  model of the given index n with r_eff ~ 5kpc at z=0.6-1.2
###
###  python py_measure-completeness.py /Users/atomczak/DATA/ORELSE/data_sg0023/noBackground_B.fits
###                                    /Users/atomczak/DATA/ORELSE/data_sg0023/Cl_0023+0423_tomczak_B_wht.fits
###                                    /Users/atomczak/DATA/ORELSE/data_sg0023/segmentation_B.fits
###                                    /Users/atomczak/DATA/ORELSE/data_sg0023/psf_B.fits
###                                    2.5
###                                    25.0
###       INPUTS:
###  (1) background subtracted image
###  (2) weight map
###  (3) segmentation map
###  (4) PSF image
###  (5) sersic index
###  (6) photometric zeropoint


import os
import sys
import mypy
import pylab
import subprocess
from astropy import wcs, modeling, cosmology, units
import scipy
from scipy import ndimage, signal
from astropy.io import fits
from photutils import CircularAperture, aperture_photometry

cosmo = cosmology.FlatLambdaCDM(Om0=0.3, H0=70)

class Sersic2D(modeling.Fittable2DModel):
    """
    Two dimensional Sersic surface brightness profile.
    Parameters
    ----------
    amplitude : float
        Central surface brightness, within r_eff.
    r_eff : float
        Effective (half-light) radius
    n : float
        Sersic Index.
    x_0 : float, optional
        x position of the center.
    y_0 : float, optional
        y position of the center.
    ellip : float, optional
        Ellipticity.
    theta : float, optional
        Rotation angle in radians, counterclockwise from
        the positive x-axis.
    See Also
    --------
    Gaussian2D, Moffat2D
    Notes
    -----
    Model formula:
    .. math::
        I(x,y) = I(r) = I_e\exp\left[-b_n\left(\frac{r}{r_{e}}\right)^{(1/n)}-1\right]
    The constant :math:`b_n` is defined such that :math:`r_e` contains half the total
    luminosity, and can be solved for numerically.
    .. math::
        \Gamma(2n) = 2\gamma (b_n,2n)
    Examples
    --------
    .. plot::
        :include-source:
        import numpy as np
        from astropy.modeling.models import Sersic2D
        import matplotlib.pyplot as plt
        x,y = pylab.meshgrid(pylab.arange(100), pylab.arange(100))
        mod = Sersic2D(amplitude = 1, r_eff = 25, n=4, x_0=50, y_0=50,
                       ellip=.5, theta=-1)
        img = mod(x, y)
        log_img = pylab.log10(img)
        plt.figure()
        plt.imshow(log_img, origin='lower', interpolation='nearest',
                   vmin=-1, vmax=2)
        plt.xlabel('x')
        plt.ylabel('y')
        cbar = plt.colorbar()
        cbar.set_label('Log Brightness', rotation=270, labelpad=25)
        cbar.set_ticks([-1, 0, 1, 2], update_ticks=True)
        plt.show()
    References
    ----------
    .. [1] http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
    """

    amplitude = modeling.Parameter(default=1)
    r_eff = modeling.Parameter(default=1)
    n = modeling.Parameter(default=4)
    x_0 = modeling.Parameter(default=0)
    y_0 = modeling.Parameter(default=0)
    ellip = modeling.Parameter(default=0)
    theta = modeling.Parameter(default=0)
    _gammaincinv = None

    @classmethod
    def evaluate(cls, x, y, amplitude, r_eff, n, x_0, y_0, ellip, theta):
        """Two dimensional Sersic profile function."""

        if cls._gammaincinv is None:
            try:
                from scipy.special import gammaincinv
                cls._gammaincinv = gammaincinv
            except ValueError:
                raise ImportError('Sersic2D model requires scipy > 0.11.')

        bn = cls._gammaincinv(2. * n, 0.5)
        a, b = r_eff, (1 - ellip) * r_eff
        cos_theta, sin_theta = pylab.cos(theta), pylab.sin(theta)
        x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
        x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
        z = pylab.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)

        return amplitude * pylab.exp(-bn * (z ** (1 / n) - 1))





###  Inputs
print '\nreading inputs...'
imname = sys.argv[1]
whtname = sys.argv[2]
segname = sys.argv[3]
psfname = sys.argv[4]
imdat = fits.getdata(imname)
imhead = fits.getheader(imname)
imwcs = wcs.WCS(imname)
whtdat = fits.getdata(whtname)
segdat = fits.getdata(segname)
psfdat = fits.getdata(psfname)
sersic_index = float(sys.argv[5])
zp = float(sys.argv[6])
print 'done!'


###  creating 2D sersic model
psf_pxscale = 0.2

xcen, ycen = pylab.array(imdat.shape) / 2
ra1, dec1 = imwcs.all_pix2world([xcen], [ycen], 1)
ra2, dec2 = imwcs.all_pix2world([xcen+1], [ycen], 1)
imdat_pxscale = mypy.radec_sep(ra1[0], dec1[0], ra2[0], dec2[0])

###  5kpc effective radius
r_eff = cosmo.arcsec_per_kpc_proper(0.9) * (5 * units.kpc)
r_eff = (r_eff / (imdat_pxscale * units.arcsec)).value

psf_width, psf_height = psfdat.shape
xgrid, ygrid = pylab.meshgrid(range(psf_width), range(psf_height))
sersic_model = Sersic2D(amplitude=1., r_eff=r_eff, n=sersic_index, 
	                    x_0=psf_width/2, y_0=psf_height/2,
	                    ellip=0., theta=0.)

sersic_model_conv = signal.convolve2d(sersic_model(xgrid, ygrid), psfdat, mode='same')
sersic_model_conv /= sersic_model_conv.sum()

fits.writeto('./psfs/sersic_n%.1f_%s' % (sersic_index, os.path.basename(imname)), sersic_model_conv)
print '\nwrote to ./psfs/sersic_n%.1f_%s' % (sersic_index, os.path.basename(imname))




###  making output directories
try: os.mkdir('./catalogs')
except: pass
try: os.mkdir('./psfs')
except: pass





###  zooming psf to pixel scale of image
psf_pxscale = 0.2

xcen, ycen = pylab.array(imdat.shape) / 2
ra1, dec1 = imwcs.all_pix2world([xcen], [ycen], 1)
ra2, dec2 = imwcs.all_pix2world([xcen+1], [ycen], 1)
imdat_pxscale = mypy.radec_sep(ra1[0], dec1[0], ra2[0], dec2[0])

psfdat_zoom = ndimage.zoom(psfdat, psf_pxscale / imdat_pxscale)

fits.writeto('./psfs/psf_%s' % os.path.basename(imname), psfdat_zoom, clobber=1)
print '\nwrote to ./psfs/psf_%s' % os.path.basename(imname)





###  estimate FWHM
dx, dy = psfdat_zoom.shape
xx, yy = pylab.meshgrid(range(-dx/2+1, dx/2+1), range(-dy/2+1, dy/2+1))
rr = (xx**2 + yy**2)**0.5

fratios = psfdat_zoom / psfdat_zoom.max()
ind = pylab.where(abs(fratios - 0.5) == abs(fratios - 0.5).min())
fwhm_px = 2 * rr[ind][0]







###  search for existing catalog of empty sky positions
try:
	sky_positions = mypy.readcat('./catalogs/sky-positions_%s.cat' % (os.path.basename(imname)[:-5]))
	xy_final = pylab.array(zip(sky_positions.x, sky_positions.y))
	noise_estimate = mypy.nmad(sky_positions.flux_aper)
except:
	print '\nidentifying empty sky positions...'

	###  finding empty sky positions, must be 3 hwhms away from nearest object
	N0 = 3 * 10**3

	xy_rand = pylab.rand(N0, 2)
	xy_rand[:,0] *= (imdat.shape[1] - 3*dx)
	xy_rand[:,1] *= (imdat.shape[0] - 3*dx)
	xy_rand[:,0] += 1.5*dx
	xy_rand[:,1] += 1.5*dx

	apers = CircularAperture(xy_rand, r=3*fwhm_px/2.)
	phot_dat = aperture_photometry(imdat, apers)
	phot_wht = aperture_photometry(whtdat, apers)
	phot_seg = aperture_photometry(segdat, apers)

	inds_poswht = pylab.find(phot_wht['aperture_sum'] > 0)
	wht = phot_wht['aperture_sum'] / pylab.median(phot_wht['aperture_sum'][inds_poswht])

	###  final sky positions with (1) positive weight and (2) not near an existing object
	inds = pylab.find((wht > 0.5) & (phot_seg['aperture_sum'] == 0))
	xy_final = xy_rand[inds]

	###  estimate aperture flux noise in positions with (1) positive weight adn (2) not near an existing object
	noise_estimate = mypy.nmad(phot_dat['aperture_sum'][inds])

	outer = open('./catalogs/sky-positions_sersic_n%.1f_%s.cat' % (sersic_index, os.path.basename(imname)[:-5]), 'w')
	outer.write('# id x y flux_aper\n')
	for ind in inds:
		outer.write('%i' % ind)
		outer.write('  %.2f' % (xy_rand[ind][0]+1))
		outer.write('  %.2f' % (xy_rand[ind][1]+1))
		outer.write('  %.4e' % phot_dat['aperture_sum'][ind])
		outer.write('\n')
	outer.close()

	print 'done!'








###  iterating magnitudes around the estimated limit
scale_factors = 2**(pylab.arange(-0.5, 5., 0.25))

magnitudes = []
frac_completeness = []

for si, scale_factor in enumerate(scale_factors):

	mock_imdat = imdat * 1.
	#mock_obj = (psfdat_zoom/psfdat_zoom.sum()) * (scale_factor*noise_estimate)
	mock_obj = (sersic_model_conv/sersic_model_conv.sum()) * (scale_factor*noise_estimate)

	mag_mock_obj = zp - 2.5 * pylab.log10(mock_obj.sum())
	magnitudes.append(mag_mock_obj)

	###  adding psf to empty sky positions
	for xi, yi in xy_final:
		xi, yi = int(xi+0.5), int(yi+0.5)
		mock_imdat[yi-dy/2:yi+(dy+1)/2, xi-dx/2:xi+(dx+1)/2] += mock_obj

	mockname = './mock_%02i_%.2f_%s' % (si+1, mag_mock_obj, os.path.basename(imname))
	fits.writeto(mockname, mock_imdat, header=imhead, clobber=1)
	print '\nwrote to %s\n' % mockname



	###  run sextractor, without assoc
	#outcat_name = './catalogs/mock_%02i_%.2f_%s_all.cat' % (si+1, mag_mock_obj, os.path.basename(imname)[:-5])
	#cmd = 'sex -c ./ddefault.sex'
	#cmd += ' %s' % mockname
	#cmd += ' -WEIGHT_IMAGE %s' % whtname
	#cmd += ' -PARAMETERS_NAME ddefault.param'
	#cmd += ' -CHECKIMAGE_TYPE NONE'
	#cmd += ' -CATALOG_NAME %s' % outcat_name
	#call = subprocess.call(cmd, shell=1)

	###  run sextractor, with assoc
	outcat_name = './catalogs/mock_%02i_%.2f_%s.cat' % (si+1, mag_mock_obj, os.path.basename(imname)[:-5])
	cmd = 'sex -c ./ddefault.sex'
	cmd += ' %s' % mockname
	cmd += ' -WEIGHT_IMAGE %s' % whtname
	cmd += ' -PARAMETERS_NAME ddefault_assoc.param'
	cmd += ' -CHECKIMAGE_TYPE NONE'
	cmd += ' -CATALOG_NAME %s' % outcat_name
	cmd += ' -ASSOC_NAME ./catalogs/sky-positions_%s.cat' % (os.path.basename(imname)[:-5])
	cmd += ' -ASSOC_PARAMS 2,3'
	cmd += ' -ASSOC_RADIUS %.1f' % fwhm_px
	call = subprocess.call(cmd, shell=1)


	###  count up the detected objects
	outcat = pylab.loadtxt(outcat_name, ndmin=2)
	n_detected = len(outcat)
	frac_completeness.append((n_detected * 1. / len(xy_final)))


	print '\n\timage:       %s' % os.path.basename(imname)[:-5]
	print '\tprogress:    %i / %i' % (si+1, len(scale_factors))
	print '\tmagnitude:   %.2f' % mag_mock_obj
	print '\tcompletness: %.2f\n' % (n_detected * 1. / len(xy_final))


	###  delete intermediate files
	os.remove(mockname)
	os.remove(outcat_name)

	###  break loop if previous 4 frac_completeness are >=0.99
	if pylab.sum(frac_completeness[-4:]) >= 4*0.99:
		break




###  writing complenetess curve to text file
outer = open('./catalogs/completeness_sersic_n%.1f_%s.cat' % (sersic_index, os.path.basename(imname)[:-5]), 'w')
outer.write('# magnitude frac_completeness\n')
for mi, fi in zip(magnitudes, frac_completeness):
	outer.write('%.3f  %.4f\n' % (mi, fi))
outer.close()





















