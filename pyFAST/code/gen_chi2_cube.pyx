
import numpy
cimport numpy

from cpython cimport array
import array


def test5(list lis):
	cdef array.array arr = array.array('d', lis)
	cdef int n  = len(arr)
	cdef int nn = 100000
	cdef double sum = 0.
	sums = []
	while nn > 0:
		nn -= 1
		while n > 0:
			n -= 1
			sum += arr[n]
		sums.append(sum)
	return sums







###  Note: For testing purposes, I'm hard-coding
###        for model_flux_cube to 3 dimensions;
###        ltau, lage, and Av
def gen_chi2_cube(numpy.ndarray flx, 
	              numpy.ndarray err, 
	              numpy.ndarray terr, 
	              numpy.ndarray model_flux_cube):

	###  reading lengths of the model cube
	cdef int n_filters = flx.shape[0]
	cdef int n_ltau = model_flux_cube.shape[0]
	cdef int n_lage = model_flux_cube.shape[1]
	cdef int n_Av   = model_flux_cube.shape[2]
	cdef int n0_filters = flx.shape[0]
	cdef int n0_ltau = model_flux_cube.shape[0]
	cdef int n0_lage = model_flux_cube.shape[1]
	cdef int n0_Av   = model_flux_cube.shape[2]


	###  creating cython arrays
	cdef array.array m
	cdef array.array fluxes = array.array('d', flx)
	cdef array.array errors = array.array('d', err)
	cdef array.array temp_errors = array.array('d', terr)

	cdef array.array etot = array.array('d', numpy.sqrt(err**2 + (terr*flx)**2))

	cdef double wmm, wfm, scale



	###  generating chi2red array
	cdef double chi2red
	scale_cube = []
	chi2red_cube = []
	while n_ltau > 0:
		n_ltau -= 1
		scale_cube.append([])
		chi2red_cube.append([])
		while n_lage > 0:
			n_lage -= 1
			scale_cube[-1].append([])
			chi2red_cube[-1].append([])
			while n_Av > 0:
				n_Av -= 1

				m = array.array('d', model_flux_cube[n0_ltau - n_ltau - 1][n0_lage - n_lage - 1][n0_Av - n_Av - 1])


				###  calculating best scale factor for model
				wmm = 0.
				wfm = 0.
				while n_filters > 0:
					n_filters -= 1
					wmm += m[n_filters] * m[n_filters] / etot[n_filters] / etot[n_filters]
					wfm += m[n_filters] * flx[n_filters] / etot[n_filters] / etot[n_filters]
				n_filters = n0_filters * 1
				scale = wfm / wmm
				scale_cube[-1][-1].append(scale)


				###  calculating chi2red of bestfit model
				chi2red = 0.
				while n_filters > 0:
					n_filters -= 1
					chi2red += (fluxes[n_filters] - m[n_filters]*scale)**2 / etot[n_filters]**2 / n0_filters
				chi2red_cube[-1][-1].append(chi2red)
				n_filters = n0_filters * 1
			n_Av = n0_Av * 1
		n_lage = n0_lage * 1
	n_ltau = n0_ltau * 1

	return scale_cube, chi2red_cube





def integrate_filter(numpy.ndarray sed_w, 
	                 numpy.ndarray sed_f,
	                 numpy.ndarray filter_w,
	                 numpy.ndarray filter_t):

	cdef array.array sed_wavelengths     = array.array('d', sed_w)
	cdef array.array sed_fluxes          = array.array('d', sed_f)
	cdef array.array filter_wavelengths  = array.array('d', filter_w)
	cdef array.array filter_transmission = array.array('d', filter_t)
	cdef int n_elements = filter_w.shape[0]

	###  interpolating sed fluxes to filter's wavelength grid
	cdef array.array sed_fluxes_interp = array.array('d', numpy.interp(filter_wavelengths, sed_wavelengths, sed_fluxes))


	###  performing trapzoidal integrations
	cdef double numerator   = 0.
	cdef double denomenator = 0.
	while n_elements > 1:
		n_elements -= 1
		mean_flux = (sed_fluxes_interp[n_elements-1] + sed_fluxes_interp[n_elements]) / 2.
		mean_trans = (filter_transmission[n_elements-1] + filter_transmission[n_elements]) / 2.
		delta_wavelength = abs(filter_wavelengths[n_elements-1] - filter_wavelengths[n_elements])

		numerator   += delta_wavelength * mean_trans * mean_flux
		denomenator += delta_wavelength * mean_trans

	return numerator / denomenator


from libc.stdlib cimport malloc, free
def integrate_filter2(numpy.ndarray sed_w, 
	                  numpy.ndarray sed_f,
	                  numpy.ndarray filter_w,
	                  numpy.ndarray filter_t):

	cdef int n1 = sed_w.shape[0]
	cdef int n2 = filter_w.shape[0]

	cdef double *sed_fluxes          = <double *>malloc(n1 * sizeof(double))
	cdef double *sed_wavelengths     = <double *>malloc(n1 * sizeof(double))
	cdef double *filter_wavelengths  = <double *>malloc(n2 * sizeof(double))
	cdef double *filter_transmission = <double *>malloc(n2 * sizeof(double))

	for i in range(n1):
		sed_fluxes[i] = sed_f[i]
		sed_wavelengths[i] = sed_w[i]

	for i in range(n2):
		filter_transmission[i] = filter_t[i]
		filter_wavelengths[i] = filter_w[i]









