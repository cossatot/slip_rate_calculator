#!/usr/bin/env python
import sys
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *

sys.path.append('../slip_rate_tools/')
import slip_rate_tools as srt

class SlipRateWindow(QMainWindow, slipRateWindow.Ui_SlipRateWindow):
    def __init__(self, parent=None):
        super(SlipRateWindow, self).__init__(parent)

        self.setupUi(self)
        self.setWindowTitle('Slip Rate Calculator')

        self.textBrowser.append("Welcome to the Slip Rate Calculator!")


        # main running buttons
        self.runButton.clicked.connect(self.run_slip_history)
        self.cancelButton.clicked.connect(self.cancel_run)
        self.plotButton.clicked.connect(self.plot_results)

        # config buttons
        self.importButton.clicked.connect(self.import_config)
        self.exportButton.clicked.connect(self.export_config)



    # main running functions
    def run_slip_history(self):
        # assemble variables from menu
        # start 
        # run (need to have some process)
        # report output (connect process to textBrowser)
        
        pass

    def cancel_run(self):
        # cancel run process

        pass

    def plot_results(self):
        # maybe plot some stuff, maybe open up a new window that
        # has lots of options before plotting
        pass


    # Config functions
    def concat_config_options(self):

    # TODO: Need to check to see how failures (empty boxes, etc.) affect
    # running of functions. Where to handle the failures? print custom errors?
        run_config = {}

        run_config['n_iters'] = float( self.nItersLineEdit.text() )
        run_config['zero_offset_age'] = float( self.zeroOffsetLineEdit.text() )
        run_config['random_seed'] = self.randSeedCheckBox.isChecked()
        run_config['random_seed_value'] = float( self.randSeedLineEdit.text() )
        run_config['force_increasing'] = self.forceIncrCheckBox.isChecked()
        run_config['slip_reversals'] = self.slipRevCheckBox.isChecked()
        run_config['fit_type'] = self.get_fit_type()
        run_config['n_linear_pieces'] = float( self.nPiecesLineEdit.text() ) 

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



app = QApplication(sys.argv)
mainWindow = SlipRateWindow()
mainWindow.show()
app.exec_()
