
import os
import array
import numpy
import subprocess
from collections import namedtuple


def readparams(name, comments='#'):

	###  gunzip if necessary
	gzipped = 0
	if name[-3:] == '.gz':
		gzipped = 1
		subprocess.call('gunzip %s' % name, shell=1)
		name = name[:-3]

	initial_read = numpy.loadtxt(name, dtype=str, comments='#')
	for i in range(len(initial_read)):
		initial_read[i][2] = initial_read[i][2].strip()
		initial_read[i][2] = initial_read[i][2].strip("'")

	custom_catalog = namedtuple('custom_catalog', initial_read[:,0])
	catalog = custom_catalog(*initial_read[:,2])

	###  re-gzip if necessary
	if gzipped:
		subprocess.call('gzip %s' % name, shell=1)

	return catalog

class res_filter:
	def __init__(self, header):
		self.name = ''
		self.err_name = ''
		self.header = header
		self.data = []
		self.lams = []
		self.trans = []
		self.l0 = 0.
def readres(res):
	'''
	#  Outputs the filter transmission curves from
	#  a *.res file commonly used with EAZY/FAST
	#
	#  res --- full path name to the *.res file
	#
	#  EXAMPLE:
	#	>>> res = '../Filters/FILTER.RES.v6.R300'
	#	>>> mypy.sps.read_filter_from_res(10, res)
	#	array([[ 2.31000e+03   0.00000e+00],
	#		   [ 2.31500e+03   0.00000e+00],
	#			  ...
	#		   [ 9.34500e+03   0.00000e+00],
	#		   [ 9.35000e+03   0.00000e+00]])
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
			filters[-1].data = numpy.array(filters[-1].data)
			filters[-1].lams = numpy.array(filters[-1].lams)
			filters[-1].trans = numpy.array(filters[-1].trans)
			filters.append(res_filter(line))

	filters[-1].data = numpy.array(filters[-1].data)
	filters[-1].lams = numpy.array(filters[-1].lams)
	filters[-1].trans = numpy.array(filters[-1].trans)
	dat.close()

	###  estimating central wavelength
	for f in filters:
		areas = []
		for i in range(len(f.lams)):
			areas.append(numpy.trapz(f.trans[:i], f.lams[:i]))
		f.l0 = numpy.interp(0.5, numpy.array(areas)/max(areas), f.lams)

	
	return filters



def read_binary(file_read, type='i', n=1, swap=False):
	line_read = array.array(type)
	line_read.fromfile(file_read, n)
	if swap: line_read.byteswap()
	if len(line_read) == 1:
		return line_read[0]
	else:
		return numpy.asarray(line_read)



def read_ised(filename):
	""" ( seds, ages, vs ) = ezgal.utils.read_ised( file )
	
	Read a bruzual and charlot binary ised file.
	
	:param file: The name of the ised file
	:type file: string
	:returns: A tuple containing model data
	:rtype: tuple
	
	.. note::
		All returned variables are numpy arrays.  ages and vs are one dimensional arrays, and seds has a shape of (vs.size,ages.size)
	
	**units**
	Returns units of:
	
	=============== ===============
	Return Variable   Units
	=============== ===============
	seds			Ergs/s/cm**2/Hz
	ages			Years
	vs			  Hz
	=============== ===============
	"""

	if not os.path.isfile(filename):
		raise ValueError('File not found: %s' % filename)

	# open the ised file
	fopen = open(filename, 'rb')

	# start reading
	junk  = read_binary(fopen)
	n_ages = read_binary(fopen)

	ages = numpy.asarray(read_binary(fopen, type='f', n=n_ages))

	# read in a bunch of stuff that I'm not interested in but which I read like this to make sure I get to the right spot in the file
	imf_lim = read_binary(fopen, type='f', n=2)
	imf_seg = read_binary(fopen, n=1)

	if imf_seg > 0:
		imf_par = numpy.zeros(imf_seg * 6)
	else:
		imf_par = 0

	for i in range(imf_seg):
		imf_par[i*6: (i+1)*6] = read_binary(fopen, type='f', n=6)

	info1  = read_binary(fopen, type='f', n=5)
	id     = read_binary(fopen, type='b', n=80)
	info2  = read_binary(fopen, type='f', n=4)
	id     = read_binary(fopen, type='b', n=80)
	id     = read_binary(fopen, type='b', n=80)
	iop    = read_binary(fopen, type='i', n=1)
	stelib = read_binary(fopen, type='i', n=1)
	junk   = read_binary(fopen)
	junk   = read_binary(fopen)
	junk   = read_binary(fopen)

	imf = [imf_lim, imf_seg, imf_par]
	info = [info1, info2, iop, stelib]
	id = str(id)

	if info1[3] == 0:
		n_extra = 12
	else:
		n_extra = 10

	# read in the wavelength data
	n_wavelengths = read_binary(fopen)

	# read wavelengths and convert to frequency (comes in as Angstroms)
	# also reverse the array so it will be sorted after converting to frequency
	wavelengths = read_binary(fopen, type='f', n=n_wavelengths)[::-1]

	# create an array for storing SED info
	seds = np.zeros((n_wavelengths, n_ages))
	ssteps = 0

	# now loop through and read in all the ages
	for i in range(n_ages):
		junk = read_binary(fopen, n=2)
		n_wavelengths = read_binary(fopen)

		seds[:,i] = read_binary(fopen, type='f', n=n_wavelengths)[::-1]
		ssteps = read_binary(fopen)
		spec = read_binary(fopen, type='f', n=ssteps)

		if i == 0:
			specs = numpy.zeros((ssteps, n_ages))
		specs[i, 0:n_wavelengths] = spec


	extras = np.zeros((n_ages, n_extra))
	for i in range(n_extra):
		junk = read_binary(fopen, n=2)
		tsteps = read_binary(fopen)
		extra = read_binary(fopen, type='f', n=n_ages)
		extras[:,i] = extra

	fopen.close()

	masses = extras[:,1]
	sfh    = extras[:,2]

	# sort by wavelength axis
	sinds = wavelengths.argsort()

	return ages, wavelengths[sinds], seds[sinds,:], ages, masses, sfh



