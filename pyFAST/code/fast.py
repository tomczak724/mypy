
import sys
import time
import read
import ezgal
import numpy

t_start = time.time()

if len(sys.argv) > 1: paramfile = sys.argv[1]
else: paramfile = 'fast.param'


params = read.readparams(paramfile)
print '\nRead parameter file: %s' % paramfile


res = read.readres(params.FILTERS_RES)
translate = numpy.loadtxt('%s.translate' % params.CATALOG, dtype=str, comments='#')
print 'Read filter files: %s, %s.translate' % (params.FILTERS_RES, params.CATALOG)





###  generating model cube/gridspace
grid_z = numpy.arange(float(params.Z_MIN), float(params.Z_MAX), float(params.Z_STEP))
grid_lage = numpy.arange(float(params.LOG_AGE_MIN), float(params.LOG_AGE_MAX), float(params.LOG_AGE_STEP))
grid_ltau = numpy.arange(float(params.LOG_TAU_MIN), float(params.LOG_TAU_MAX), float(params.LOG_TAU_STEP))
grid_Av = numpy.arange(float(params.A_V_MIN), float(params.A_V_MAX), float(params.A_V_STEP))






path2models = '%s/ised_%s.%s' % (params.LIBRARY_DIR, params.SFH, params.RESOLUTION)


model_cube = []

for i_ltau in range(len(grid_ltau)):
	ltau = grid_ltau[i_ltau]
	modelname = '%s/%s_%s_%s_z%s_ltau%.1f.ised' % (path2models, params.LIBRARY, params.RESOLUTION, params.IMF, params.METAL[2:], ltau)
	model_cube.append(ezgal.model(modelname))














