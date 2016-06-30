
import glob
import time
import numpy
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects


###  Enter your name
inspector_name = raw_input('\n\tEnter you last name: ')
output_file_name = 'visual_inspection_%s.txt' % inspector_name
output_file = open(output_file_name, 'w')



###  Grabbing mugshot files
mugshots = glob.glob('./data/mugshot*png')



###  Initializing interactive figure
fig = pyplot.figure(figsize=(9.1875, 8.875))
subplot = fig.add_subplot(111, label='subplot')
fig.subplots_adjust(left=0., top=1., right=1., bottom=0.)
subplot.xaxis.set_visible(0)
subplot.yaxis.set_visible(0)
subplot.axhline(0, color='gray', lw=1.5)
subplot.axis([-0.5, 1088.5, -180, 890.5])


###  Spacing paramters for buttons
gap = 0.06
button_width = (1. - 5*gap) / 4.
button_height = 0.06
bottom_offset = 0.04

buttons = []
button_names = ['all_good', 'missed_star', 'bad_fit', 'other']
button_colors = ['lime', 'orange', 'red', 'gray']

output_file.write('# id')
for bname in button_names:
	output_file.write(' %s' % bname)


###  Adding buttions to plot window
for bi in range(len(button_names)):
	buttons.append(subplot.figure.add_axes([(bi+1)*gap+bi*button_width, bottom_offset, button_width, button_height], axisbg='w', label=button_names[bi]))
	buttons[-1].patch.set_facecolor(button_colors[bi])
	buttons[-1].patch.set_alpha(0.6)
	buttons[-1].xaxis.set_visible(0)
	buttons[-1].yaxis.set_visible(0)

	t = buttons[-1].text(0.5, 0.5, button_names[bi], transform=buttons[-1].transAxes, fontsize=18, fontweight='bold', verticalalignment='center', horizontalalignment='center')



###  Plotting first mugshot
i_gal = 0
image_data = pyplot.imread(mugshots[i_gal])
image_graphics = subplot.imshow(image_data[::-1])

progress_text = subplot.text(0.03, bottom_offset + 1.3*button_height, 'Progress:  1 / %i' % len(mugshots), 
	                         transform=subplot.transAxes, color='gray', fontsize=16,
	                         path_effects=[PathEffects.withStroke(linewidth=2., foreground='k')])



###  This script gets executed every time you click somewhere in the plot window
def onclick(event):
	global i_gal
	button_flags = numpy.zeros(len(buttons), dtype=int)


	###  grab the name of which subplot was clicked
	clicked_subplot = event.inaxes.get_label()


	###  if the subplot is one of the buttons, do this
	if clicked_subplot in button_names:
		
		button_flags[button_names.index(clicked_subplot)] += 1

		field = mugshots[i_gal].split('_')[1]
		galaxy_id = int(mugshots[i_gal].split('_')[0].split('mugshot')[1])

		###  write the parameters to the output file
		output_file.write('\n%7i ' % galaxy_id)
		for flag in button_flags:
			output_file.write(' %i' % flag)

		###  increment by +1 and plot the next mugshot, or quit if no more left
		i_gal += 1
		if i_gal >= len(mugshots):
			subplot.text(0.5, 0.5, '\n  Thanks for your time!                               \n', transform=subplot.transAxes, fontsize=35, fontweight='bold', verticalalignment='center', horizontalalignment='center', bbox={'facecolor':'#ffff80', 'linewidth':3})
			pyplot.draw()
			time.sleep(4)

			pyplot.close()
			output_file.write('\n')
			output_file.close()
			print '\n\n\twrote to: %s\n' % output_file_name
		else:
			redraw(i_gal)


###  script to plot the next mugshot
def redraw(i):
	progress_text.set_text('Progress:  %i / %i' % (i+1, len(mugshots)))
	image_data = pyplot.imread(mugshots[i])
	image_graphics.set_data(image_data[::-1])
	pyplot.draw()


###  IMPORTANT --- this initiates a connection between the plot window and the onclick() script
cid = fig.canvas.mpl_connect('button_press_event', onclick)
pyplot.show()










