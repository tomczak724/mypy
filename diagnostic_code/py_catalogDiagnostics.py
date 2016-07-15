
import mypy
import pylab
import subprocess


def readcat_gzip(catname):
	subprocess.call('gunzip %s.gz' % catname, shell=1)
	output = mypy.readcat(catname)
	subprocess.call('gzip %s' % catname, shell=1)
	return output

def readzout_gzip(catname):
	subprocess.call('gunzip %s.gz' % catname, shell=1)
	output = mypy.readzout(catname)
	subprocess.call('gzip %s' % catname, shell=1)
	return output

def plot_uvj_region(subplot, z, color='r', lw=2, ls='-', zorder=99):
    
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

    p = subplot.plot( (-10,c1), (floor,floor), color=color, lw=lw, ls=ls, zorder=zorder )
    p = subplot.plot( (c1,wall), (floor,c2), color=color, lw=lw, ls=ls, zorder=zorder )
    p = subplot.plot( (wall,wall), (c2,10), color=color, lw=lw, ls=ls, zorder=zorder )




print '\nreading catalogs...'
'''
version = 'cl0910+5422_v0.0.2'
cat_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs/cl0910+5422_v0.0.2'
cat = readcat_gzip('%s/%s.cat' % (cat_dir, version))
zout = readzout_gzip('%s/%s.zout' % (cat_dir, version))
fout = readcat_gzip('%s/%s.fout' % (cat_dir, version))
rf = readcat_gzip('%s/%s.restframe' % (cat_dir, version))
fiducial_filter = 0
filternames = ['V', 'Rc', 'I+', 'Z+', 'J', 'K']
fluxes = [[cat.fluxaper_V, cat.fluxauto_V],
          [cat.fluxaper_Rc, cat.fluxauto_Rc],
          [cat.fluxaper_Iplus, cat.fluxauto_Iplus],
          [cat.fluxaper_Zplus, cat.fluxauto_Zplus],
          [cat.fluxaper_J, cat.fluxauto_J],
          [cat.fluxaper_K, cat.fluxauto_K]]
'''

version = 'rcs0224-0002_v0.0.2'
cat_dir = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/RCS_0224/release/rcs0224-0002_v0.0.2'
cat = readcat_gzip('%s/%s.cat' % (cat_dir, version))
zout = readzout_gzip('%s/%s.zout' % (cat_dir, version))
fout = readcat_gzip('%s/%s.fout' % (cat_dir, version))
rf = readcat_gzip('%s/%s.restframe' % (cat_dir, version))
fiducial_filter = 1
filternames = ['B', 'V', 'R+', 'I+', 'Z+', 'J', 'K']
fluxes = [[cat.fluxaper_B, cat.fluxauto_B],
          [cat.fluxaper_V, cat.fluxauto_V],
          [cat.fluxaper_Rplus, cat.fluxauto_Rplus],
          [cat.fluxaper_Iplus, cat.fluxauto_Iplus],
          [cat.fluxaper_Zplus, cat.fluxauto_Zplus],
          [cat.fluxaper_J, cat.fluxauto_J],
          [cat.fluxaper_K, cat.fluxauto_K]]

print 'done!\n'





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
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)]) + 0.2
	asdf = sp1.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=pylab.cm.gray_r)

	pylab.savefig('colorDist%02i_%s-%s_%s.png' % (i+1, filternames[fiducial_filter], filtername, version))
	sp1.clear()

pylab.close()









###################################
###  UVJ and (U-V) vs M* diagrams
###################################

zbins = [[0.50, 0.75],
         [0.75, 1.00],
         [1.00, 1.50]]

fig = pylab.figure(figsize=(12.3, 13.7))
sp2 = fig.add_subplot(222)
sp1 = fig.add_subplot(221)
sp4 = fig.add_subplot(224)
sp3 = fig.add_subplot(223)
sps = [sp1, sp2, sp3, sp4]
fig.subplots_adjust(wspace=0, hspace=0, top=0.96)

for i in range(3):

	zlo, zhi = zbins[i]

	for sp in sps:
		sp.minorticks_on()
		sp.grid()

	sp1.set_title('%.2f < z < %.2f' % (zlo, zhi), size=22)
	sp2.set_title('%s' % version, size=22)

	sp3.set_xlabel('( V - J )$_{rest}$', size=24)
	sp1.set_ylabel('( U - V )$_{rest}$', size=24)
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=24)
	sp3.set_ylabel('( U - V )$_{rest}$', size=24)

	sp1.axis([-0.7, 2.7, -0.7, 2.7])
	sp2.axis([7.3, 11.8, -0.7, 2.7])
	sp3.axis([-0.7, 2.7, -0.7, 2.7])
	sp4.axis([7.3, 11.8, -0.7, 2.7])

	plot_uvj_region(sp1, (zlo+zhi)/2.)
	plot_uvj_region(sp3, (zlo+zhi)/2.)




	### photometric sample
	inds = pylab.find((cat.use == 1) &
		              (fout.z > zlo) & 
		              (fout.z < zhi) &
		              (zout.odds > 0.8) &
		              (fout.lmass > 0))

	uv = -2.5 * pylab.log10(rf.restflux_U[inds] / rf.restflux_V[inds])
	vj = -2.5 * pylab.log10(rf.restflux_V[inds] / rf.restflux_J[inds])
	lmass = fout.lmass[inds]

	###  UVJ diagram
	hist2d, xedges, yedges = pylab.histogram2d(vj, uv, bins=(100, 100), range=([-0.7, 2.7], [-0.7, 2.7]))
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)]) + 0.2
	asdf = sp1.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=pylab.cm.gray_r)
	sp1.set_aspect('auto')

	###  (U-V) vs M* diagram
	hist2d, xedges, yedges = pylab.histogram2d(lmass, uv, bins=(100, 100), range=([7.3, 11.8], [-0.7, 2.7]))
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)]) + 0.2
	asdf = sp2.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=pylab.cm.gray_r)
	sp2.set_aspect('auto')

	t = sp3.text(1, 1.05, 'photometric:  N = %i ' % len(inds), transform=sp3.transAxes, verticalalignment='center', horizontalalignment='center', bbox={'facecolor':'yellow'})




	### FULL spectroscopic sample
	inds = pylab.find((cat.use == 1) &
		              (cat.z_spec > 0) &
		              (fout.lmass > 0))

	uv = -2.5 * pylab.log10(rf.restflux_U[inds] / rf.restflux_V[inds])
	vj = -2.5 * pylab.log10(rf.restflux_V[inds] / rf.restflux_J[inds])
	lmass = fout.lmass[inds]

	sp3.plot(vj, uv, 'ko', ms=1)
	sp4.plot(lmass, uv, 'ko', ms=1)

	t = sp3.text(1, 0.95, 'spectroscopic:  N = %i ' % len(inds), transform=sp3.transAxes, verticalalignment='center', horizontalalignment='center', bbox={'facecolor':'yellow'})


	pylab.savefig('UVJmass_%.2fz%.2f_%s.png' % (zlo, zhi, version))
	sp1.clear()
	sp2.clear()
	sp3.clear()
	sp4.clear()

pylab.close()














