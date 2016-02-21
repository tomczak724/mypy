import os
import pylab as pl
import numpy as np
import pyfits as pf
import subprocess as sp
import matplotlib as mpl
from scipy import interpolate
from scipy import optimize as opt

c = 2.99792458e18  # Ang / s


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
    if splineinterp > 1:
        lamax = np.linspace( filter[0][0], filter[-1][0], int(len(filter)*splineinterp) )
        trans = interpolate.spline( filter[:,0], filter[:,1], lamax )
        filter = np.array( zip(lamax,trans) )
        
    filter_normed = filter * 1.
    filter_normed[:,1] /= np.trapz( filter[:,1], filter[:,0] )

    fconvolved = filter_normed[:,1] * np.interp( filter_normed[:,0], spec[:,0], spec[:,1] )

    return np.trapz( fconvolved, filter_normed[:,0] )




def uvj_select(uv, vj, z):
    """
    #  INPUT:
    #    uv  --->  Restframe (U-V) color (in AB system)  [array]
    #    vj  --->  Restframe (V-J) color (in AB system)  [array]
    #    z  ---->  Galaxy redshifts  [array]
    #
    #  OUTPUT: quiescent=0, star-forming=1
    #
    #  EXAMPLE:
    #            [ SF, SF, Qui, SF,  ....  Qui, Qui]
    #     array( [ 1,  1,  0,   1,   ....  0,   0 ] )
    #
    #  based on Whitaker+2011
    """
    if type(uv)==list or type(uv)==float or type(uv)==int or type(uv)==long or type(uv)==np.float64: uv = np.array(uv)
    if type(vj)==list or type(vj)==float or type(vj)==int or type(vj)==long or type(vj)==np.float64: vj = np.array(vj)

    floor = np.zeros(len(z))
    wall = np.zeros(len(z))
    slope = np.zeros(len(z))
    intercept = np.zeros(len(z))

    floor[ np.where( (0.0<=z) & (z<1.5) ) ] = 1.3
    floor[ np.where( (1.5<=z) & (z<2.0) ) ] = 1.3
    floor[ np.where( (2.0<=z) ) ] = 1.2

    wall[ np.where( (0.0<=z) & (z<1.5) ) ] = 1.6
    wall[ np.where( (1.5<=z) & (z<2.0) ) ] = 1.5
    wall[ np.where( (2.0<=z) ) ] = 1.4

    slope[ np.where( z<0.5 ) ] = 0.88
    slope[ np.where( 0.5<=z ) ] = 0.88

    intercept[ np.where( z<0.5 ) ] = 0.69
    intercept[ np.where( 0.5<=z ) ] = 0.59

    outer = np.zeros(len(z))
    outer[ np.where( (uv<slope*vj+intercept) | (uv<floor) | (vj>wall) )[0] ] = 1
    return outer



def uvj_select_region( z, plot=True, kolor='k', lw=3, ls='-', zorder=99):
    """
    #  DESCRIPTION:
    #    Returns the defined region for selecting
    #    star-forming vs quiescent galaxies from
    #    Whitaker+2011. Will plot it unless told
    #    not to.
    #
    #  INPUT:
    #    z  --->  desired redshift
    #
    #  OUTPUT:
    #    floor  ----->  minimum U-V color to be quiescent
    #    corner1  --->  V-J coordinate at turning point at floor
    #    corner2  --->  U-V coordinate at turning point at wall
    #    wall  ------>  maximum V-J color to be quiescent
    """
    
    if 0<=z<1.5:
        floor = 1.3
        wall = 1.6
    elif 1.5<=z<2.0:
        floor = 1.3
        wall = 1.5
    elif 2.0<=z:
        floor = 1.2
        wall = 1.4

    if 0<z<0.5:
        slope = 0.88
        intercept = 0.69
    elif 0.5<=z:
        slope = 0.88
        intercept = 0.59

    c1 = (floor - intercept) / slope
    c2 = slope*wall + intercept

    if plot:
        p = pl.plot( (-10,c1), (floor,floor), color=kolor, lw=lw, ls=ls, zorder=zorder )
        p = pl.plot( (c1,wall), (floor,c2), color=kolor, lw=lw, ls=ls, zorder=zorder )
        p = pl.plot( (wall,wall), (c2,10), color=kolor, lw=lw, ls=ls, zorder=zorder )

    return floor, c1, c2, wall




class feature():
   def __init__(self, lam, name, type):
      self.lam = lam      # approx wavelength
      self.name = name    # name of the feature
      self.type = type    # absorbtion, emission or sky

def spec_lines( lam, z=0 ):
    """
  #  This script prints to the screen various
  #  spectral features nearby the specified
  #  wavelength. If supplied a redshift, it
  #  will also print the redshifted wavelength.
  #
  #  An input of lam=0 will print all lines.
    """
    lines = [ feature( float(i[0]), i[1], i[2] ) for i in np.loadtxt('/Users/atomczak/.spectral_lines',dtype=str) ]

    s = "\n\n\tRest-lam(Ang)\t"
    if z>0: s += "lam at z="+str(z)+"\t"
    s += " Name\t      Type\n"
    print(s)

    found = False
    for line in lines:

        if lam==0: 
            s = "\t  "+ str(line.lam)[:7] +"\t"
            if z>0: s += "  "+str(line.lam*(1+z))[:7] +"\t"
            s += " "+line.name +"\t      "+ line.type
            print(s)
        elif abs((lam-line.lam)/lam)<0.05:
            found = True

            s = "\t  "+ str(line.lam)[:7] +"\t"
            if z>0: s += "  "+str(line.lam*(1+z))[:7] +"\t"
            s += " "+line.name +"\t      "+ line.type
            print(s)

    print("")
    if not found: print("\t\t  SORRY... NO MATCHES FOUND\n")



def centroid( image, x, y, fwhm, sky=0 ):
    """
    #
    #  DESCRIPTION:
    #    Computes centroid(s) by fitting gaussians in X,Y
    #
    #  INPUTS:
    #    image = 2D array of image data
    #    x,y   = nominal coordinates of star(s)
    #    fwhm  = approximate full-width-half-max
    #    sky   = background sky value to be subtracted.
    #            if not specified, it is assumed to be 0.
    #
    #  OUTPUTS:
    #    xout,yout = calulated centroid coordinates
    #
    #  NOTES:
    #    - outputs (and inputs) are treated as indices.
    #      that is... x,y will be offset by 1 than in ds9.
    #
    #  Author: A.R. Tomczak
    #
    """
    if type(x)==int or type(x)==float or type(x)==np.float64: x = [x]
    if type(y)==int or type(y)==float or type(y)==np.float64: y = [y]

###  Checking sizes of inputs
    if len(x)!=len(y): raise IOError("Input X and Y do not match")

    xout,yout = [],[]
    coords = [ [round(x[i]),round(y[i])] for i in range(len(x)) ]

###  Looping over stars
    for xyi in range(len(x)):
        x0i, y0i = x[xyi], y[xyi]
        xtmp, ytmp = coords[xyi]

#xxx
#        print x0i, y0i

###  Setting # iterations to find peak pixel to 2*fwhm
###  If it takes longer, script will give up.
        found_peak_pix = False
        for tolerance in range(int(round(2*fwhm))):
            cutout = image[ytmp-1:ytmp+2, \
                           xtmp-1:xtmp+2]
            for i in range(3):
                for j in range(3):
                    if cutout[j][i] == cutout.max(): break
                if cutout[j][i] == cutout.max(): break
            if i==1 and j==1:
                found_peak_pix = True
                break
            else:
                xtmp += (i-1)
                ytmp += (j-1)
        if not found_peak_pix: raise IOError("Could not find peak pixel")

        
###  Now to fit gaussians
        xax = np.arange( int(xtmp-fwhm*2/3),int(xtmp+fwhm*2/3)+1 )
        yax = np.arange( int(ytmp-fwhm*2/3),int(ytmp+fwhm*2/3)+1 )
        xflux = image[ ytmp, xax ]
        yflux = image[ yax, xtmp ]


        xguess = [max(xflux) - sky, xtmp, fwhm/4., sky]
        yguess = [max(yflux) - sky, ytmp, fwhm/4., sky]

        xsigma = (abs(xax - xtmp) + 1)**0.5
        ysigma = (abs(xax - xtmp) + 1)**0.5

        xfit, xcov = opt.curve_fit(mypy.gauss1d, xax, xflux, p0=xguess, sigma=xsigma)
        yfit, ycov = opt.curve_fit(mypy.gauss1d, yax, yflux, p0=yguess, sigma=ysigma)
        xout.append(xfit[1])
        yout.append(yfit[1])

###  Done! Now to spit it out.
    if len(xout)==1: return xout[0], yout[0]
    else: return np.array(zip(xout, yout))



def progress(n, ntot, percent=10):
    """
    #  If n is an integer percentage in units of
    #  "percent" then will print to the screen.
    #
    #  EXAMPLES:
    #
    #  >>> progress(2, 8, 25)
    #  >>>   PROGRESS: 25 %
    #  >>>
    #  >>> progress(3, 8, 25)
    #  >>> progress(4, 8, 25)
    #  >>>   PROGRESS: 50 %
    """
    ratio = n * 100. / ntot
    prog_array = np.arange(int(ntot * percent / 100.), ntot+1, int(ntot * percent / 100.))
    if int(n) in prog_array:
        return (1, " PROGRESS: " + str(int(ratio+0.5)) + " %")
    else:
        return (0, " PROGRESS: " + str(int(ratio+0.5)) + " %")




def get_Hz(wavelength):
    '''
    #  DESCRIPTION
    #    Returns the frequency for the given
    #    wavelength (in Angstroms)
    '''
    return c / wavelength

def get_lambda(frequency):
    '''
    #  DESCRIPTION
    #    Returns the wavelength (in Angstrom)
    #     for the given frequency.
    '''
    return c / frequency





class res_filter:
    def __init__(self, header):
        self.header = header
        self.data = []
        self.lams = []
        self.trans = []


def read_res(res):
    '''
    #  Outputs the filter transmission curves from
    #  a *.res file commonly used with EAZY/FAST
    #
    #  res --- full path name to the *.res file
    #
    #  EXAMPLE:
    #    >>> res = '../Filters/FILTER.RES.v6.R300'
    #    >>> mypy.sps.read_filter_from_res(10, res)
    #    array([[ 2.31000e+03   0.00000e+00],
    #           [ 2.31500e+03   0.00000e+00],
    #              ...
    #           [ 9.34500e+03   0.00000e+00],
    #           [ 9.35000e+03   0.00000e+00]])
    '''

    dat = open(res, 'r')
    lines = dat.readlines()

    filters = []

    line0 = lines[0]
    n = int(line0.split()[0])
    filters.append(res_filter(line0))

    cnt = 0
    for line in lines[1:]:
        
        if cnt < n:
            cnt, lam, tran = line.split()
            cnt = int(cnt)
            lam = float(lam)
            tran = float(tran)
            filters[-1].data.append([lam, tran])
            filters[-1].lams.append(lam)
            filters[-1].trans.append(tran)

        else:
            cnt = 0
            n = int(line.split()[0])
            filters[-1].data = np.array(filters[-1].data)
            filters[-1].lams = np.array(filters[-1].lams)
            filters[-1].trans = np.array(filters[-1].trans)
            filters.append(res_filter(line))

    filters[-1].data = np.array(filters[-1].data)
    filters[-1].lams = np.array(filters[-1].lams)
    filters[-1].trans = np.array(filters[-1].trans)
    dat.close()
    
    return filters



res = read_res('%s/data/filters.res' % os.getcwd())
v_band = res[141]

vega = np.loadtxt('%s/data/vega.dat' % os.getcwd())

ab_nu = np.array(zip(vega[:,0], np.zeros(len(vega)) + synphot(vega, v_band.data)))

def vega2ab(filt):
    '''
    #  Returns the (AB-Vega) for a given filter.
    #
    #  filt = np.array([lam_1, trans_1],
    #                   lam_2, trans_2],
    #                      ...
    #                   lam_n, trans_n]])
    '''
    f_vega = synphot(vega_nu, filt)
    f_ab = synphot(ab_nu, filt)
    return 2.5 * np.log10(f_ab / f_vega)

