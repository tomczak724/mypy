
import numpy

def lum2flux(z, cosmology, ZP_AB=25.):

	zp_c_fl = 10**((48.57 + ZP_AB) / 2.5)
	dist    = cosmology.luminosity_distance(z).value
	l2f     = 10**(numpy.log10(3.826e1) + numpy.log10(zp_c_fl) - numpy.log10(4. * numpy.pi * (1. + z) * (dist*3.0856776e8)**2))

	return l2f

