#!/usr/bin/env python
import sys
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *

try:
    from qtconsole.rich_ipython_widget import RichIPythonWidget
    from qtconsole.inprocess import QtInProcessKernelManager
except ImportError:
    from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
    from IPython.qt.inprocess import QtInProcessKernelManager

import slipRateWindow

sys.path.append('../slip_rate_tools/')
import slip_rate_tools as srt

import pandas as pd
import numpy as np
from collections import OrderedDict
import json

import qt_plots as qtp

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

        # main running buttons
        self.runButton.clicked.connect(self.run_slip_history)
        self.cancelButton.clicked.connect(self.cancel_run)
        self.plotButton.clicked.connect(self.plot_results)

        # config buttons
        self.importButton.clicked.connect(self.import_config)
        self.exportButton.clicked.connect(self.export_config)

        # IPython console
        self.console = EmbedIPython(srt=srt, plt=plt, qtp=qtp)
        self.vlayout_for_ipython.addWidget(self.console)


        # offset marker table
        self.tabledata = test_table_data

        self.tablemodel = OffsetMarkerTableModel(self.tabledata, 
                                                 offset_table_header)

        self.offsetMarkerTableView.setModel(self.tablemodel)
        self.offsetMarkerTableView.horizontalHeader()

        self.addOffsetMarkerButton.clicked.connect(self.add_offset_marker_row)
        self.removeOffsetMarkerButton.clicked.connect(
                                                self.remove_offset_marker_row)

        self.offsetMarkerTableView.model().layoutChanged.connect(
                                                         self.push_table_data)
        self.push({'table_header': offset_table_header})
        

    def push(self, obj):
        '''convenience function for pushing objects to the console'''
        self.console.kernel.shell.push(obj)

    def add_offset_marker_row(self):
        self.tabledata.append(
            ['name', 0., 'age_type', 0., 'age_err_type', 'age_unit', 
                     0., 'offset_type', 0., 'offset_err_type', 'offset_unit'])

        self.offsetMarkerTableView.model().layoutChanged.emit()

    def remove_offset_marker_row(self):
        selections = self.offsetMarkerTableView.selectedIndexes()[0]
        row_index = selections.row()
        rm = self.offsetMarkerTableView.model().removeRow(row_index)

        self.offsetMarkerTableView.model().layoutChanged.emit()

    def push_table_data(self):
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
                                        #'offset_list':offset_list,
                                        'tabledata':self.tabledata})
        
        self.console.execute(
            'offset_list = srt.offset_list_from_gui(tabledata, table_header)')

        self.console.execute(
            'res_df, age_arr, offset_arr, n_pieces_best = '
            + 'srt.run_interp_from_gui(offset_list, rc)')

        #self.console.kernel.shell.push({'off_table':self.tabledata})

        #pass

    def cancel_run(self):
        # cancel run process
        # none of the options seem to work...

        #self.console.execute("\x03\r\n")
        #self.console.interrupt_kernel()
        self.console.kernel_client.stop_channels()
        # pass

    def plot_results(self):

        #plot_cmd = ('qtp.results_plots_for_gui(res_df, age_arr, rc, '
        plot_cmd = ('plot_win = qtp.results_plots_for_gui(res_df, age_arr, rc,'
                                                   +' n_pieces_best,'
                                                   +' offset_arr, offset_list,'
                                                   +' show_samples=True)' )
        self.console.execute(plot_cmd)

        return


    # Config functions
    def concat_config_options(self):

    # TODO: Need to check to see how failures (empty boxes, etc.) affect
    # running of functions. Where to handle the failures? print custom errors?
        run_config = OrderedDict()

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

    def concat_all_variables(self):
        all_vars = OrderedDict()

        all_vars['run_config'] = self.concat_config_options()
        all_vars['offset_markers'] = self.tabledata_to_dict()

        return all_vars

    def dict_to_config(self, input_dict):
        if 'run_config' in input_dict.keys():
            dic = input_dict['run_config']
        else:
            dic = input_dict

        self.nItersLineEdit.setText( str(dic['n_iters']))
        self.zeroOffsetLineEdit.setText( str(dic['zero_offset_age']))
        self.randSeedCheckBox.setChecked( dic['random_seed'])
        self.randSeedLineEdit.setText( str(dic['random_seed_value']))
        self.forceIncrCheckBox.setChecked( dic['force_increasing'])
        self.slipRevCheckBox.setChecked( dic['slip_reversals'])
        if dic['fit_type'] == 'linear':
            self.linearFitRadio.setChecked(True)
        elif dic['fit_type'] == 'piecewise':
            self.piecewiseFitRadio.setChecked(True)
        elif dic['fit_type'] == 'cubic':
            self.cubicFitRadio.setChecked(True)
        #self.nPiecesSpin
        return

    def import_config(self):
        infile = QFileDialog.getOpenFileName(self, filter='*.json')

        fp = open(infile, 'r')
        
        all_vars = json.load(fp, object_pairs_hook=OrderedDict)
        
        self.dict_to_tabledata(all_vars)
        self.dict_to_config(all_vars)
        
        return

    def export_config(self):
        
        all_vars = self.concat_all_variables()

        output_name = QFileDialog.getSaveFileName(self, 
                                                  caption='Save File As',
                                                  filter='*.json')
        with open(output_name, 'w+') as f:
            json.dump(all_vars, f, indent=2)
        
        return

    def tabledata_to_dict(self):

        offset_marker_d = OrderedDict()

        for row in self.tabledata:
            name = row[0]

            offset_marker_d[name] = srt.offset_marker_dict_from_row(row[1:],
                                                    offset_table_header[1:])
        return offset_marker_d

    def dict_to_tabledata(self, input_dict):

        offset_list = []
       
        if 'offset_markers' in input_dict.keys():
            dic = input_dict['offset_markers']

        else:
            dic = input_dict

        for key in dic:
            row_list = [key]

            for head in offset_table_header[1:]:
                row_list.append(dic[key][head])

            offset_list.append(row_list)
        
        self.tabledata = offset_list
        
        self.tablemodel = OffsetMarkerTableModel(self.tabledata, 
                                                 offset_table_header)
        
        self.offsetMarkerTableView.setModel(self.tablemodel)

        return





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


    def removeRow(self, row, parent=QModelIndex()):
        self.removeRows(row, 1, parent)

    def removeRows(self, row, count, parent=QModelIndex()):
        self.beginRemoveRows(parent, row, row+count-1)
        for i in reversed(range(count)):
            self.arraydata.pop(row+i)
        self.endRemoveRows()
        return True




offset_table_header = ['Name', 'Age', 'Age_Type', 'Age_Err', 'Age_Err_Type',
                       'Age_Units', 'Offset', 'Offset_Type', 'Offset_Err', 
                       'Offset_Err_Type', 'Offset_Units']


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
test_table_data = [['T1', 24., 'mean', 8., 'sd', 'ka', 
                    list(t1.offset_m), 'list', list(t1.rel_prob), 'probs', 'm'],
                   ['Qa', 50., 'mean', 20., 'sd', 'ka', 
                    75., 'mean', 20., 'sd', 'm'],
                   ['Qao', 100., 'mean', 32., 'sd', 'ka',
                    132., 'mean', 33., 'sd', 'm']]



app = QApplication(sys.argv)
mainWindow = SlipRateWindow()
mainWindow.show()
app.exec_()
