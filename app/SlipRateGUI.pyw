#!/usr/bin/env python
import sys
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
from IPython.qt.inprocess import QtInProcessKernelManager


import slipRateWindow

sys.path.append('../slip_rate_tools/')
import slip_rate_tools as srt

import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


'''
Test data.  This will be removed when the adding data interface is built.
'''

offset_df = pd.read_csv('../test_data/offsets.csv')
offset_df['offset_m'] = offset_df.offset_in * 200.

t1 = offset_df[offset_df.unit == 'T1']
qa = offset_df[offset_df.unit == 'Qa']
qao = offset_df[offset_df.unit == 'Qao']

#qa['offset_m'] += 200.

t1_age = {'mean': 24., 'sd':8.}
qa_age = {'mean': 50., 'sd':20.}
qao_age = {'mean':100., 'sd':32.}

qao_age['mean'] += 200

T1 = srt.OffsetMarker(age_mean=t1_age['mean'], age_sd=t1_age['sd'],
                      offset_vals=t1.offset_m, offset_probs=t1.rel_prob)

Qa = srt.OffsetMarker(age_mean=qa_age['mean'], age_sd=qa_age['sd'],
                      offset_vals=qa.offset_m, offset_probs=qa.rel_prob)

Qao = srt.OffsetMarker(age_mean=qao_age['mean'], age_sd=qao_age['sd'],
                      offset_vals=qao.offset_m, offset_probs=qao.rel_prob)


offset_list = [T1, Qa, Qao]

plt=''



class EmbedIPython(RichIPythonWidget):

    def __init__(self, **kwarg):
        super(RichIPythonWidget, self).__init__()
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel()
        self.kernel = self.kernel_manager.kernel
        self.kernel.gui = 'qt4'
        self.kernel.shell.push(kwarg)
        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()



class SlipRateWindow(QMainWindow, slipRateWindow.Ui_MainWindow):
    def __init__(self, parent=None):
        super(SlipRateWindow, self).__init__(parent)

        self.setupUi(self)
        self.setWindowTitle('Slip Rate Calculator')

        #self.textBrowser.append("Welcome to the Slip Rate Calculator!")


        # main running buttons
        self.runButton.clicked.connect(self.run_slip_history)
        self.cancelButton.clicked.connect(self.cancel_run)
        self.plotButton.clicked.connect(self.plot_results)

        # config buttons
        self.importButton.clicked.connect(self.import_config)
        self.exportButton.clicked.connect(self.export_config)

        # IPython console
        self.console = EmbedIPython(srt=srt, plt=plt)
        self.vlayout_for_ipython.addWidget(self.console)


        # offset marker table
        self.tabledata = test_table_data

        self.tablemodel = OffsetMarkerTableModel(self.tabledata, 
                                                 offset_table_header)
        self.offsetMarkerTableView.setModel(self.tablemodel)
        self.offsetMarkerTableView.horizontalHeader()

        self.addOffsetMarkerButton.clicked.connect(self.add_offset_marker_row)
        

    def push(self, **kwarg):
        '''convenience function for pushing objects to the console'''
        self.console.kernel.shell.push(kwarg)

    def add_offset_marker_row(self):
        self.tabledata.append(['name', 0., 'age_type', 0., 'age_err_type'])
        self.offsetMarkerTableView.model().layoutChanged.emit()

        print(self.tabledata)
        print(type(self.tabledata))
        self.push({'tabledata':self.tabledata})




    # main running functions
    def run_slip_history(self):
        # assemble variables from menu
        # start 
        # run (need to have some process)
        # report output (connect process to textBrowser)
        #self.console.execute("pwd")
        
        run_config_dict = self.concat_config_options()
        self.console.kernel.shell.push({'rc':run_config_dict, 
                                        'offset_list':offset_list})

        self.console.execute(
            'res_df, age_arr, offset_arr = srt.run_interp_from_gui('
            + 'offset_list, rc)')

        #self.console.kernel.shell.push({'off_table':self.tabledata})

        #pass

    def cancel_run(self):
        # cancel run process

        #self.console.execute("\x03\r\n")
        #self.console.interrupt_kernel()
        self.console.kernel_client.stop_channels()
        # pass

    def plot_results(self):
        # maybe plot some stuff, maybe open up a new window that
        # has lots of options before plotting

        #plot_cmd = 'srt.plot_histograms_from_gui(rc, res_df)'
        plot_cmd = ('srt.plot_slip_histories_from_gui(res_df, age_arr, rc,'
                                                    +'offset_arr, offset_list,'
                                                    +'show_samples=True)' )
                                                     

        self.console.execute(plot_cmd)

        #self.console.kernel.shell.p

        #plt.hist(np.random.random(7000), bins=50)
        #plt.show()

        #pass


    # Config functions
    def concat_config_options(self):

    # TODO: Need to check to see how failures (empty boxes, etc.) affect
    # running of functions. Where to handle the failures? print custom errors?
        run_config = {}

        run_config['n_iters'] = int( self.nItersLineEdit.text() )
        run_config['zero_offset_age'] = float( self.zeroOffsetLineEdit.text() )
        run_config['random_seed'] = self.randSeedCheckBox.isChecked()
        run_config['random_seed_value']=int(float(self.randSeedLineEdit.text()))
        run_config['force_increasing'] = self.forceIncrCheckBox.isChecked()
        run_config['slip_reversals'] = self.slipRevCheckBox.isChecked()
        run_config['fit_type'] = self.get_fit_type()
        run_config['n_linear_pieces'] = int( self.nPiecesSpinBox.value() ) 

        return run_config

    def get_fit_type(self):

        if self.linearFitRadio.isChecked():
            fit_type = 'linear'
        elif self.piecewiseFitRadio.isChecked():
            fit_type = 'piecewise'
        elif self.cubicFitRadio.isChecked():
            fit_type = 'cubic'

        return fit_type

    def concat_offset_markers(self):
        # don't really know how to deal with this right now
        pass

    def concat_all_variables(self):
        all_vars = {}

        for key, val in concat_config_options():
            all_vars[key] = val

        for key, val in concat_offset_markers():
            all_vars[key] = val

        return all_vars

    def import_config(self):
        # open file browser, select file
        # try reading as csv, json
        # output run_config dict
        pass

    def export_config(self):
        # open file browser
        # select either csv or json
        # write to appropriate file, don't allow others
        
        all_vars = concat_all_variables()
        
        pass


class OffsetMarkerTableModel(QAbstractTableModel):
    '''docs'''
    def __init__(self, data_in, header_in, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.arraydata = data_in
        self.header = header_in


    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != Qt.DisplayRole:
            return None
        return self.arraydata[index.row()][index.column()] 

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        if len(self.arraydata) > 0:
            return len(self.arraydata[0])
        else:
            return 0

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.header[col]

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            self.arraydata[index.row()][index.column()] = value
            return True
        else:
            return False

    def flags(self, index):
        return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable










offset_table_header = ['Name', 'Age', 'Age_Type', 'Age_Err', 'Age_Err_Type']
                      # 'Age_Units', 'Offset', 'Offset_Type', 'Offset_Err', 
                      # 'Offset_Err_Type', 'Offset_Units']


test_table_data = [['T1', 24., 'mean', 8., 'sd'],
                   ['Qa', 50., 'mean', 20., 'sd'],
                   ['Qao', 100., 'mean', 32., 'sd']]



app = QApplication(sys.argv)
mainWindow = SlipRateWindow()
mainWindow.show()
app.exec_()
