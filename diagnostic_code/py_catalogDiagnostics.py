
import mypy
import pylab
import subprocess


def readcat_gzip(catname):
	subprocess.call('gunzip %s.gz' % catname, shell=1)
	output = mypy.readcat(catname)
	subprocess.call('gzip %s' % catname, shell=1)
	return output



version = 'cl0910+5422_v0.0.2'
cat_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs/cl0910+5422_v0.0.2'

print '\nreading catalogs...'
cat = readcat_gzip('%s/%s.cat' % (cat_dir, version))
fout = readcat_gzip('%s/%s.fout' % (cat_dir, version))
rf = readcat_gzip('%s/%s.restframe' % (cat_dir, version))
print 'done!\n'

fiducial_filter = 0
filternames = ['V', 'Rc', 'I+', 'Z+', 'J', 'K']
fluxes = [[cat.fluxaper_V, cat.fluxauto_V],
          [cat.fluxaper_Rc, cat.fluxauto_Rc],
          [cat.fluxaper_Iplus, cat.fluxauto_Iplus],
          [cat.fluxaper_Zplus, cat.fluxauto_Zplus],
          [cat.fluxaper_J, cat.fluxauto_J],
          [cat.fluxaper_K, cat.fluxauto_K]]




########################################
###  Observed color magnitude diagrams
########################################

fig = pylab.figure()
sp1 = fig.add_subplot(111)
fig.subplots_adjust(top=0.96)

for i, filtername in enumerate(filternames):
	if filtername == filternames[fiducial_filter]:
		continue

	sp1.set_title(version)
	sp1.minorticks_on()
	axlims = [18, 28, -2, 5]
	sp1.axis(axlims)
	sp1.grid()
	sp1.set_xlabel('%s$_{AUTO}$' % filternames[fiducial_filter])
	sp1.set_ylabel('( %s - %s )$_{APER}$' % (filternames[fiducial_filter], filtername))

	inds = pylab.find((cat.use == 1) & 
		              (fluxes[i][0] > 0) &
		              (fluxes[i][1] > 0) &
		              (fluxes[fiducial_filter][0] > 0) &
		              (fluxes[fiducial_filter][1] > 0))


	color = -2.5 * pylab.log10(fluxes[fiducial_filter][0][inds] / fluxes[i][0][inds])
	magnitude = 25 - 2.5 * pylab.log10(fluxes[fiducial_filter][1][inds])

	hist2d, xedges, yedges = pylab.histogram2d(magnitude, color, bins=(500, 500))
	extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]]
	hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])
	asdf = sp1.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=pylab.cm.gray_r)

	pylab.savefig('colorDist%02i_%s-%s_%s.png' % (i+1, filternames[fiducial_filter], filtername, version))
	sp1.clear()

pylab.close()









########################################
###  Observed color magnitude diagrams
########################################

zbins = [[0.50, 0.75],
         [0.75, 1.00],
         [1.00, 1.50]]

fig = pylab.figure(figsize=(8.1, 8.8))
sp2 = fig.add_subplot(322)
sp1 = fig.add_subplot(321)
sp4 = fig.add_subplot(324)
sp3 = fig.add_subplot(323)
sp6 = fig.add_subplot(326)
sp5 = fig.add_subplot(325)
sps = [[sp1, sp2], [sp3, sp4], [sp5, sp6]]
fig.subplots_adjust(wspace=0, hspace=0)

for i in range(3):

	zlo, zhi = zbins[i]

	spa = sps[i][0]
	spb = sps[i][1]

	spa.minorticks_on()
	spb.minorticks_on()
	spa.grid()
	spb.grid()

	spa.axis([-0.7, 2.7, -0.7, 2.7])
	spb.axis([7.2, 11.8, -0.7, 2.7])

	spa.set_xlabel('( V - J )$_{rest}$')
	spa.set_ylabel('( U - V )$_{rest}$')
	spb.set_xlabel('log( M* / M$_{\odot}$ )')
	spb.set_ylabel('( U - V )$_{rest}$')

	spa.text(0.03, 0.87, '%.2f < z < %.2f' % (zlo, zhi), transform=spa.transAxes)
	spb.text(0.03, 0.87, '%.2f < z < %.2f' % (zlo, zhi), transform=spb.transAxes)


	inds = pylab.find((cat.use == 1) &
		              (fout.z > zlo) & 
		              (fout.z < zhi) &
		              (fout.lmass > 0))


	uv = -2.5 * pylab.log10(rf.restflux_U[inds] / rf.restflux_V[inds])
	vj = -2.5 * pylab.log10(rf.restflux_V[inds] / rf.restflux_J[inds])
	lmass = fout.lmass[inds]

	spa.plot(vj, uv, 'ko', ms=1)
	spb.plot(lmass, uv, 'ko', ms=1)





pylab.savefig('colorDist%02i_%s-%s_%s.png' % (i+1, filternames[fiducial_filter], filtername, version))














