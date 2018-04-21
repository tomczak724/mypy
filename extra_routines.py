
import sys
import time
import numpy


def progress_bar(n, ntot, print_start=True):
    """
    Prints a basic progress bar to show progress of
    time-consuming "for" loops

    Parameters
    ----------
    n : int
        Iteration number of loop

    ntot : int
        Total number of interations

    print_start : bool
        Print date+time of iteration 0
    """

    ###  print date+time of first iteration                                 
    if print_start and n == 0:
        print('\n   start: %s' % time.ctime())

    prog = (n + 1) * 100. / ntot
    text = '\r   progress  %2i %%  [ ' % prog

    ###  each "=" corresponds to 2% of progress                             
    n_units = int(prog / 2)
    text += '=' * n_units
    text += ' ' * (50 - n_units)
    text += ' ]'

    sys.stdout.write(text)
    sys.stdout.flush()


def mad(arr):
    """
    Returns the Median Absolute Deviation
    """
    abs_dev = abs(numpy.array(arr) - numpy.median(arr))
    return number.median(abs_dev)


def nmad(arr):
    """
    Returns the Normalized Median Absolute Deviation
    """
    return 1.483 * mad(arr)
