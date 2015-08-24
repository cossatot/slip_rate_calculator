import numpy as np
import pandas as pd
#from sklearn.neighbors import KernelDensity
from slip_rate_tools import *

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib import collections as mc
from matplotlib import gridspec

#from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from PyQt4.QtGui import *
from PyQt4.QtCore import *



def results_plots_for_gui(res_df, age_arr, run_config_dict, n_pieces,
                          offset_arr=None, offset_list=None,
                          show_data=False, show_samples=True):

 
    # get appropriate symbols for number of iterations
    n_iters = run_config_dict['n_iters']

    if n_iters < 10:
        sym = 'o'
    elif 10 <= n_iters < 100:
        sym = '.'
    elif 100 <= n_iters:
        sym = ','

    nbins = (int(np.log10(n_iters))+1) * 10
    #canvas = MplCanvas(num_subplots=2, num_pieces=n_pieces)
    canvas = PlotWindow(num_subplots=2, num_pieces=n_pieces)

    slip_history_ax = canvas.ax1
    slip_history_ax.set_ylabel('Modern offset (m)')

    # rate and rate change histograms
    if n_pieces == 1:
        histograms_ax = canvas.ax3
        histograms_ax.hist(res_df.m, bins=nbins)
        histograms_ax.set_xlabel('Rate (mm/yr)')
        #histograms_ax.set_ylabel('count')
        slip_history_ax.set_xlabel('Age of geologic feature')

    elif n_pieces == 2:
        slip_history_ax.tick_params( axis='x', labelbottom='off')
    
        rate_histories = rate_history_fits(res_df, age_arr, 
                                           run_config_dict, n_pieces)
        ymax = res_df[['m1', 'm2']].values.max()
        ymin = res_df[['m1', 'm2']].values.min()
        rate_history_ax = canvas.ax2
        rate_history_ax.set_xlabel('Time (ka)')
        rate_history_ax.set_ylabel('Rate (mm/yr)')
        rate_history_ax.add_collection(rate_histories)
        rate_history_ax.set_ylim([ymin, ymax])

        m1_ax = canvas.ax3
        m1_ax.hist(res_df.m1, bins=nbins)
        #m1_ax.set_xlabel('Younger slip rate (mm/yr)')
        m1_ax.tick_params( axis='x', labelbottom='off')

        bkpt_ax = canvas.ax4a
        bkpt_ax.hist(res_df.breakpt, bins=nbins)
        bkpt_ax.set_xlabel('Timing of rate change (ka)')

        rate_change_ax = canvas.ax4b
        rate_change_ax.hist((res_df.m1 - res_df.m2), bins=nbins)
        rate_change_ax.set_xlabel('Rate change (mm/yr)')

        scatter_ax = canvas.ax5
        scatter_ax.plot(res_df.m1, res_df.m2, sym)
        scatter_ax.set_xlabel('Younger slip rate (mm/yr)')
        scatter_ax.set_ylabel('Older slip rate (mm/yr)')
        
        m2_ax = canvas.ax6
        m2_ax.hist(res_df.m2, orientation='horizontal', bins=nbins)
        #m2_ax.yaxis.set_visible(False)
        m2_ax.tick_params( axis='y', labelbottom='off')


    if show_data == True:
        #ax.errorbar() # need to finish
        pass


    if show_samples == True:

        # TODO: loop through offset markers for different colors here
        #slip_history_ax.plot(age_arr.ravel(), offset_arr.ravel(), sym)
        
        for col in range(age_arr.shape[1])[1:]:
            slip_history_ax.plot(age_arr[:,col], offset_arr[:,col], sym)

    slip_histories = slip_history_fits(res_df, age_arr, run_config_dict, 
                                       n_pieces)

    slip_history_ax.add_collection(slip_histories)
    slip_history_ax.set_xlim(left=0.)
    slip_history_ax.set_ylim(bottom=0.)
    slip_history_ax.invert_xaxis()

    
    canvas.show()


def line_thickness_adjust(n_lines, exp=0.4, numer=1):
    return numer / n_lines**exp


def slip_history_fits(res_df, age_arr, run_config_dict, n_pieces):

    n_iters = age_arr.shape[0]
    lw = line_thickness_adjust(n_iters)
    
    if run_config_dict['fit_type'] in ('linear', 'piecewise'):
        pass
    
    elif run_config_dict['fit_type'] == 'cubic':
        raise Exception('cubic spline plotting not implemented yet')

    else:
        raise Exception('fit type needs to be linear, piecewise, or cubic')

    line_pts = get_history_line_pts_from_results(res_df, age_arr, n_pieces)

    line_coll = mc.LineCollection(line_pts, linewidths=0.05,
                                  #alpha=0.5, 
                                  colors='grey')
    return line_coll


def rate_history_fits(res_df, age_arr, run_config_dict, n_pieces):

    n_iters = age_arr.shape[0]
    lw = line_thickness_adjust(n_iters)
    

    if run_config_dict['fit_type'] in ('linear', 'piecewise'):
        pass
    
    elif run_config_dict['fit_type'] == 'cubic':
        raise Exception('cubic spline plotting not implemented yet')

    else:
        raise Exception('fit type needs to be linear, piecewise, or cubic')

    line_pts = get_rate_line_pts_from_results(res_df, age_arr, n_pieces)

    line_coll = mc.LineCollection(line_pts, linewidths=0.5,
                                  alpha=lw, colors='b')

    return line_coll


def get_rate_line_pts_from_results(res_df, age_array, n_pieces):

    if age_array.shape[0] != res_df.shape[0]:
        age_array = trim_age_offset_arrays(res_df, age_array)
    
    n_rows = res_df.shape[0]
    n_pts = n_pieces * 2

    x_array = np.zeros((n_rows, n_pts))
    y_array = np.zeros((n_rows, n_pts))
    
    if n_pieces == 1:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = age_array.max(axis=1)

        y_array[:,0] = res_df.m
        y_array[:,1] = res_df.m

    elif n_pieces == 2:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = res_df.breakpt.values
        x_array[:,2] = res_df.breakpt.values
        x_array[:,3] = age_array.max(axis=1)

        y_array[:,0] = res_df.m1
        y_array[:,1] = res_df.m1
        y_array[:,2] = res_df.m2
        y_array[:,3] = res_df.m2
    
    line_pts = [tuple(zip(x_array[i,:], y_array[i,:])) for i in range(n_rows)]

    return line_pts


def get_history_line_pts_from_results(res_df, age_array, n_pieces):
    
    if age_array.shape[0] != res_df.shape[0]:
        age_array = trim_age_offset_arrays(res_df, age_array)

    n_pts = n_pieces + 1
    n_rows = res_df.shape[0]

    x_array = np.zeros((n_rows, n_pts))
    y_array = x_array.copy()
    
    if n_pieces == 1:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = age_array.max(axis=1)

        y_array[:,0] = 0.
        y_array[:,1] = age_array.max(axis=1) * res_df.m

    elif n_pieces == 2:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = res_df.breakpt.values
        x_array[:,2] = age_array.max(axis=1)

        y_array[:,0] = 0.
        y_array[:,1] = res_df.breakpt.values * res_df.m1
        y_array[:,2] = cumulative_offsets(x_array[:,1], res_df.m1, 
                                          x_array[:,2], res_df.m2)
    else:
        raise Exception('only 1 or 2 piece lines currently implemented.')

    line_pts = [tuple(zip(x_array[i,:], y_array[i,:])) for i in range(n_rows)]

    return line_pts


class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=12, height=6, 
                 num_subplots=2, num_pieces=1):#, dpi=100):
        #super(MplCanvas, self).__init__(self)

        self.fig = Figure(figsize=(width, height),
                     #, dpi=dpi
                     )

        super(MplCanvas, self).__init__(self.fig)#, width=width, height=height, 
                                        #num_subplots=num_subplots, 
                                        #num_pieces=num_pieces)
        self.canvas = FigureCanvas(self.fig)

        if num_subplots == 1:
            self.axes = self.fig.add_subplot(111)

        elif num_subplots == 2:

            if num_pieces == 1:
                gs = gridspec.GridSpec(1,3)

                ax_1 = gs.new_subplotspec((0, 0), colspan=2)
                self.ax1 = self.fig.add_subplot(ax_1)
                
                ax_3 = gs.new_subplotspec((0, 2))
                self.ax3 = self.fig.add_subplot(ax_3)

            elif num_pieces == 2:
                gr = 4
                gc = 4
                gs = gridspec.GridSpec(gr,gc)
                    
                ax_1 = gs.new_subplotspec((0,0), colspan=gc//2, rowspan=gr//2)
                ax_2 = gs.new_subplotspec((gr//2,0), colspan=gc//2, 
                                                     rowspan=gr//2)

                self.ax1 = self.fig.add_subplot(ax_1)
                self.ax2 = self.fig.add_subplot(ax_2, sharex=self.ax1)
    
                # m1 (youngest)
                ax_3 = gs.new_subplotspec((0,gc//2), colspan=gc//4, 
                                                     rowspan=gr//2)
                self.ax3 = self.fig.add_subplot(ax_3)

                # breakpt_hist
                ax_4a = gs.new_subplotspec((0,gc-gc//4), colspan=gc//4,
                                                         rowspan=gr//4)
                self.ax4a = self.fig.add_subplot(ax_4a)
                
                # rate change hist
                ax_4b = gs.new_subplotspec((1,gc-gc//4), colspan=gc//4,
                                                         rowspan=gr//4)
                self.ax4b = self.fig.add_subplot(ax_4b)

                # m1, m2 scatter
                ax_5 = gs.new_subplotspec((gr//2, gc//2), colspan=gc//4,
                                                          rowspan=gr//2)
                self.ax5 = self.fig.add_subplot(ax_5, sharex=self.ax3)

                # m2 (older) hist
                ax_6 = gs.new_subplotspec((gr//2, gc-gc//4), colspan=gc//4,
                                                             rowspan=gr//2)
                self.ax6 = self.fig.add_subplot(ax_6, sharey=self.ax5)
        
                self.fig.subplots_adjust(hspace=0.45, bottom=0.1,
                                         left=0.05, right=0.95)


class PlotWindow(QDialog, MplCanvas):
    def __init__(self, parent=None, width=14, height=7, num_subplots=2,
                 num_pieces=1):
        super(PlotWindow, self).__init__(parent, width=width, height=height, 
                                         num_subplots=num_subplots, 
                                         num_pieces=num_pieces)

        self.main_frame = QWidget()
        self.canvas.setParent(self.main_frame)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame,
                                             coordinates=False)

        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.canvas)
        self.vbox.addWidget(self.mpl_toolbar)
        self.setLayout(self.vbox)

        

