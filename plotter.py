
###  My custom routines related to figures and plotting

import numpy


def restore_logaxes_labels(subplot, xaxis=True, yaxis=True):
    '''
    This script reformats the default labeling scheme
    for logarithmic axis labels to the "regular" format.
    e.g. axis labels of "10**-1" and "10**3" will be
    changed to "0.1" and "1000" respectively.
    '''
    xticks = subplot.get_xticks() 
    yticks = subplot.get_yticks() 

    if xaxis:
        xticks_new = []
for xi in xticks:
    base = numpy.log10(xi)
    if base >= 0:
        xi_new = '%i' % xi
    else:
        formatter = '%.' + str(int(abs(base))) + 'f'
        xi_new = formatter % xi
    xticks_new.append(xi_new)
subplot.xaxis.set_ticklabels(xticks_new)

if yaxis:
    yticks_new = []
for yi in yticks:
    base = numpy.log10(yi)
    if base >= 0:
        yi_new = '%i' % yi
    else:
        formatter = '%.' + str(int(abs(base))) + 'f'
        yi_new = formatter % yi
    yticks_new.append(yi_new)
subplot.yaxis.set_ticklabels(yticks_new)

return subplot




def add_inset(subplot, rect=[0.4, 0.4, 0.2, 0.2]):
    '''
    This script creates an Axes instance within the
    coordinate frame of the provided subplot.
    '''

    ###  coordinates of the subplot in the figure window's coordinate frame
    box = subplot.get_position()
    xlo, xhi = box.x0, box.x1
    ylo, yhi = box.y0, box.y1
    dx, dy = (xhi - xlo), (yhi - ylo)

    ###  width/height of requested axes in the figure window's coordinate frame
    sub_dx = dx * rect[2]
    sub_dy = dy * rect[3]

    ###  position of requested axes in the figure window's coordinate frame
    sub_xlo = xlo + dx * rect[0]
    sub_ylo = ylo + dy * rect[1]

    inset = subplot.figure.add_axes([sub_xlo, sub_ylo, sub_dx, sub_dy])
    return inset

