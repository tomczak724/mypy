import numpy

def calzetti(wavelength):
	R_v   = 4.05

	k_lam = numpy.zeros(len(wavelength))
	x     = wavelength / (1.e4)  # wavelength in microns
	x1    = numpy.where(x <  0.63)[0]
	x2    = numpy.where(x >= 0.63)[0]

	k_lam[x1] = 2.659 * (-2.156 + 1.509/x[x1] - 0.198/(x[x1]**2) + 0.011/(x[x1]**3)) + R_v
	k_lam[x2] = 2.659 * (-1.857 + 1.040/x[x2]) + R_v

	return k_lam / R_v

