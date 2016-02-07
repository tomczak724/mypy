
import mypy
import pylab as pl


res = mypy.sps.read_res('filters.res')


pl.savetxt('./data/HST/WFC/F435W.txt', res[1-1].data)
pl.savetxt('./data/HST/WFC/F475W.txt', res[2-1].data)
pl.savetxt('./data/HST/WFC/F555W.txt', res[3-1].data)
pl.savetxt('./data/HST/WFC/F606W.txt', res[4-1].data)
pl.savetxt('./data/HST/WFC/F775W.txt', res[5-1].data)
pl.savetxt('./data/HST/WFC/F814W.txt', res[6-1].data)
pl.savetxt('./data/HST/WFC/F850LP.txt', res[7-1].data)

pl.savetxt('./data/HST/WFPC2/F300W.txt', res[10-1].data)
pl.savetxt('./data/HST/WFPC2/F336W.txt', res[11-1].data)
pl.savetxt('./data/HST/WFPC2/F450W.txt', res[12-1].data)
pl.savetxt('./data/HST/WFPC2/F555W.txt', res[13-1].data)
pl.savetxt('./data/HST/WFPC2/F606W.txt', res[14-1].data)
pl.savetxt('./data/HST/WFPC2/F702W.txt', res[15-1].data)
pl.savetxt('./data/HST/WFPC2/F814W.txt', res[16-1].data)
pl.savetxt('./data/HST/WFPC2/F850LP.txt', res[17-1].data)

pl.savetxt('./data/HST/WFC3/F098M.txt', res[201-1].data)
pl.savetxt('./data/HST/WFC3/F105W.txt', res[202-1].data)
pl.savetxt('./data/HST/WFC3/F125W.txt', res[203-1].data)
pl.savetxt('./data/HST/WFC3/F140W.txt', res[204-1].data)
pl.savetxt('./data/HST/WFC3/F160W.txt', res[205-1].data)

pl.savetxt('./data/Spitzer/IRAC/3.6um.txt', res[18-1].data)
pl.savetxt('./data/Spitzer/IRAC/4.5um.txt', res[19-1].data)
pl.savetxt('./data/Spitzer/IRAC/5.8um.txt', res[20-1].data)
pl.savetxt('./data/Spitzer/IRAC/8.0um.txt', res[21-1].data)

pl.savetxt('./data/SDSS/u.txt', res[156-1].data)
pl.savetxt('./data/SDSS/g.txt', res[157-1].data)
pl.savetxt('./data/SDSS/r.txt', res[158-1].data)
pl.savetxt('./data/SDSS/i.txt', res[159-1].data)
pl.savetxt('./data/SDSS/z.txt', res[160-1].data)

pl.savetxt('./data/2MASS/J.txt', res[161-1].data)
pl.savetxt('./data/2MASS/H.txt', res[162-1].data)
pl.savetxt('./data/2MASS/Ks.txt', res[163-1].data)

pl.savetxt('./data/FOURSTAR/J.txt', res[248-1].data)
pl.savetxt('./data/FOURSTAR/J1.txt', res[249-1].data)
pl.savetxt('./data/FOURSTAR/J2.txt', res[250-1].data)
pl.savetxt('./data/FOURSTAR/J3.txt', res[251-1].data)
pl.savetxt('./data/FOURSTAR/H.txt', res[251-1].data)
pl.savetxt('./data/FOURSTAR/Hs.txt', res[254-1].data)
pl.savetxt('./data/FOURSTAR/Hl.txt', res[253-1].data)
pl.savetxt('./data/FOURSTAR/Ks.txt', res[254-1].data)

pl.savetxt('./data/CFHT/MegaCam/u.txt', res[88-1].data)
pl.savetxt('./data/CFHT/MegaCam/g.txt', res[89-1].data)
pl.savetxt('./data/CFHT/MegaCam/r.txt', res[90-1].data)
pl.savetxt('./data/CFHT/MegaCam/i.txt', res[91-1].data)
pl.savetxt('./data/CFHT/MegaCam/z.txt', res[92-1].data)



















