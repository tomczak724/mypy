
import numpy
from scipy import interpolate

def get_filter_central_wavelength(wavelengths, transmission):
    areas = []
    for i in range(len(wavelengths)):
        areas.append(numpy.trapz(transmission[:i], wavelengths[:i]))
    return numpy.interp(0.5, numpy.array(areas) / max(areas), wavelengths)


def synphot(spec, filter, splineinterp=1):
    """
    #  Basic script to measure integrated flux from
    #  a given spectrum in a given bandpass/filter.
    #  Flux must be in F_lambda units.
    #
    #  spec   = array( [ [lam_1,flux_1], [lam_2,flux_2], ... [lamn,flux_n] ] )
    #  filter = array( [ [lam_1,tran_1], [lam_2,tran_2], ... [lamn,tran_n] ] )
    #
    #  If splineinterp is >1 then the filter transmission
    #  curve will be subsampled by the factor splineinterp
    #  using a spline interpolation.
    """

    if splineinterp>1:
        lamax = numpy.linspace(filter[0][0], filter[-1][0], int(len(filter)*splineinterp))
        trans = interpolate.spline( filter[:,0], filter[:,1], lamax)
        filter = numpy.array(zip(lamax,trans))
        
    filter_normed = filter
    filter_normed[:,1] /= numpy.trapz(filter[:,1], filter[:,0])

    fconvolved = filter_normed[:,1] * numpy.interp(filter_normed[:,0], spec[:,0], spec[:,1])

    return numpy.trapz(fconvolved, filter_normed[:,0])



