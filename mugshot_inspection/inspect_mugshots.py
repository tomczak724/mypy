
import glob
import time
import numpy
from matplotlib import pyplot


inspector_name = raw_input('\n\tEnter you last name: ')
output_file_name = 'visual_inspection_%s.txt' % inspector_name
output_file = open(output_file_name, 'w')


mugshots = glob.glob('./data/mugshot*png')





fig = pyplot.figure(figsize=(11.0625, 10.775))
subplot = fig.add_subplot(111, label='subplot')
fig.subplots_adjust(left=0., top=1., right=1., bottom=0.)
subplot.xaxis.set_visible(0)
subplot.yaxis.set_visible(0)

###  prepping the canvas
subplot.axis([-0.5, 1088.5, -150, 890.5])


gap = 0.06
button_width = (1. - 5*gap) / 4.
button_height = 0.06
bottom_offset = 0.06


buttons = []
button_names = ['all_good', 'missed_star', 'bad_fit', 'other']
button_colors = ['lime', 'orange', 'red', 'gray']

output_file.write('# id')
for bname in button_names:
	output_file.write(' %s' % bname)

for bi in range(len(button_names)):
	buttons.append(subplot.figure.add_axes([(bi+1)*gap+bi*button_width, bottom_offset, button_width, button_height], axisbg='w', label=button_names[bi]))
	buttons[-1].patch.set_facecolor(button_colors[bi])
	buttons[-1].patch.set_alpha(0.6)
	buttons[-1].xaxis.set_visible(0)
	buttons[-1].yaxis.set_visible(0)

	t = buttons[-1].text(0.5, 0.5, button_names[bi], transform=buttons[-1].transAxes, fontsize=18, fontweight='bold', verticalalignment='center', horizontalalignment='center')




image_data = pyplot.imread(mugshots[0])
subplot.imshow(image_data[::-1])




i_gal = 0
def onclick(event):
	global i_gal

	clicked_subplot = event.inaxes.get_label()

	button_flags = numpy.zeros(len(buttons), dtype=int)
	if clicked_subplot in button_names:
		
		button_flags[button_names.index(clicked_subplot)] += 1

		field = mugshots[i_gal].split('_')[1]
		galaxy_id = int(mugshots[i_gal].split('_')[0].split('mugshot')[1])

		output_file.write('\n%7i ' % galaxy_id)
		for flag in button_flags:
			output_file.write(' %i' % flag)





		i_gal += 1
		if i_gal >= len(mugshots):
			subplot.text(0.5, 0.5, '\n  Thanks for your time!    \n', transform=subplot.transAxes, fontsize=35, fontweight='bold', verticalalignment='center', horizontalalignment='center', bbox={'facecolor':'w', 'linewidth':3})
			pyplot.draw()
			time.sleep(4)

			pyplot.close()
			output_file.write('\n')
			output_file.close()
			print '\n\n\twrote to: %s\n' % output_file_name
		else:
			redraw(i_gal)


def redraw(i):
	image_data = pyplot.imread(mugshots[i])
	subplot.imshow(image_data[::-1])
	pyplot.draw()




cid = fig.canvas.mpl_connect('button_press_event', onclick)



pyplot.show()





















'''

class clickabel_window:

	def __init__(self):
		
		self.fig = pyplot.figure(figsize=(11.0625, 10.775))
		self.subplot = self.fig.add_subplot(111)
		self.fig.subplots_adjust(left=0., top=1., right=1., bottom=0.15)
		self.subplot.xaxis.set_visible(0)
		self.subplot.yaxis.set_visible(0)

		self.fig.canvas.mpl_connect('button_press_event', self.onclick)


	def display_image_data(self, image_display, image_data):
		image_display.set_data(image_data)


	def onclick(self, event):
		x = event.xdata
		y = event.ydata
		print x, y




mugshots = glob.glob('./data/mugshot*png')


gui = clickabel_window()


for i, mugshot in enumerate(mugshots):

	image_data = pyplot.imread(mugshot)
	if i == 0:
		image_display = gui.subplot.imshow(image_data[::-1])
	else:
		gui.display_image_data(image_display, image_data[::-1])

	break
	time.sleep(3)


'''





















































