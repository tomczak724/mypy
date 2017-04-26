
import sys
import mypy
import numpy
from astropy.io import fits
from matplotlib import pyplot
pyplot.ion()


image_name = sys.argv[1]






image = fits.open(image_name)
naxis1 = image[0].header['NAXIS1']
naxis2 = image[0].header['NAXIS2']









class image_display(object):

	def __init__(self, subplot, data):

		self.subplot = subplot
		self.data = data

		###  setting defaults
		self.x_center = data.shape[1] / 2
		self.y_center = data.shape[0] / 2

		box_size1 = data.shape[1] / 10
		box_size2 = data.shape[0] / 10
		self.box_size = max([box_size1, box_size2])

		initial_cutout = data[self.y_center-self.box_size/2:self.y_center+self.box_size/2,
		                      self.x_center-self.box_size/2:self.x_center+self.box_size/2]
		bg_median = numpy.median(initial_cutout)
		bg_nmad = mypy.nmad(initial_cutout)
		self.vmin = bg_median-3*bg_nmad
		self.vmax = bg_median+5*bg_nmad

		self.display = pyplot.imshow(initial_cutout, cmap=pyplot.cm.gray,
			                         vmin=self.vmin, vmax=self.vmax)
		cid_onClick    = fig.canvas.mpl_connect('button_press_event', onClick)
		cid_onScroll   = fig.canvas.mpl_connect('scroll_event',       onScroll)
		cid_onKeyPress = fig.canvas.mpl_connect('key_press_event',    onKeyPress)



	def redraw(self):
		cutout = self.data[self.y_center-self.box_size/2:self.y_center+self.box_size/2,
		                   self.x_center-self.box_size/2:self.x_center+self.box_size/2]
		self.display.set_data(cutout)
		self.display.set_clim(self.vmin, self.vmax)
		self.display.set_extent([-0.5, self.box_size-0.5, -0.5, self.box_size-0.5])
		pyplot.draw()





def onClick(event):

	###  try to grab the name of which subplot was clicked
	try:
		chosen_subplot = event.inaxes.get_label()
	except:
		chosen_subplot = 'out_of_bounds'


	###  if the subplot is one of the buttons, do this
	if chosen_subplot == 'image_subplot':

		x_offset = int(image_subplot.box_size / 2 - event.xdata)
		y_offset = int(image_subplot.box_size / 2 - event.ydata)
		image_subplot.x_center -= x_offset
		image_subplot.y_center -= y_offset
		image_subplot.redraw()



def onScroll(event):

	###  try to grab the name of which subplot was clicked
	try:
		chosen_subplot = event.inaxes.get_label()
	except:
		chosen_subplot = 'out_of_bounds'


	###  if the subplot is one of the buttons, do this
	if chosen_subplot == 'image_subplot':

		if event.step < 0:
			fzoom = 1.05
		else:
			fzoom = 0.95

		image_subplot.box_size = int(image_subplot.box_size * fzoom)
		image_subplot.redraw()


def onKeyPress(event):

	###  try to grab the name of which subplot was clicked
	try:
		chosen_subplot = event.inaxes.get_label()
	except:
		chosen_subplot = 'out_of_bounds'


	###  if the subplot is one of the buttons, do this
	if chosen_subplot == 'image_subplot':

		###  print image stats near key press
		if event.key == 's':

			dx = 5

			x_offset = int(image_subplot.box_size / 2 - event.xdata)
			y_offset = int(image_subplot.box_size / 2 - event.ydata)

			cutout_x_center = image_subplot.x_center - x_offset
			cutout_y_center = image_subplot.y_center - y_offset

			cutout = image_subplot.data[cutout_y_center-dx:cutout_y_center+dx,
			                            cutout_x_center-dx:cutout_x_center+dx]

			mean = numpy.mean(cutout)
			median = numpy.median(cutout)
			stdev = numpy.std(cutout)
			nmad = mypy.nmad(cutout)

			###  print x, y, dx, mean, std, median, nmad
			output  = '  %6.1f' % cutout_x_center
			output += '  %6.1f' % cutout_y_center
			output += '  %2i' % (2*dx)
			output += '  %.4e' % mean
			output += '  %.4e' % stdev
			output += '  %.4e' % median
			output += '  %.4e' % nmad
			print output







###  Initializing interactive figure
fig = pyplot.figure()
fig.set_facecolor('w')

image_subplot = image_display(fig.add_subplot(111, label='image_subplot'),
	                          image[0].data)

































