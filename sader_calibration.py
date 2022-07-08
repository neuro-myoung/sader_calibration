## Import packages required to run program
from PyQt5 import QtWidgets, QtCore
from buildv3 import Ui_MainWindow
import sys

import nidaqmx
import numpy as np
from nidaqmx.constants import TerminalConfiguration
from nidaqmx.constants import VoltageUnits
from time import sleep

from matplotlib.widgets import SpanSelector

from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings
from scipy import signal

from sympy import besselk
from mpmath import *


## Generate main window class
class mywindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(mywindow, self).__init__()
        
        ## Call ui from buildui.py 
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.screen = QtWidgets.QDesktopWidget().screenGeometry()
        self.setGeometry(QtCore.QRect(50, 50, 1200,800))


        ## Link appropriate functions to widget signal events 
        self.ui.acquirebutton.clicked.connect(self.acquire)
        self.ui.selectButton.clicked.connect(self.data_select)
        self.ui.fitButton.clicked.connect(self.fit_data)
        self.ui.calculatebutton.clicked.connect(self.calculate_stiffness)
        self.ui.actionSave_Raw_Data.setShortcut("Ctrl+S")
        self.ui.actionSave_Raw_Data.triggered.connect(self.save_raw_data)
        self.ui.actionSave_Power_Spectra.triggered.connect(self.save_Pv_data)

        
    ###################### Functions #########################################
    def acquire(self):
        """ Begin acquisition of analog input 0 signal from NIDAQ device.
        This function will acquire the raw input coming through the ai0 channel
        of Dev1 using the nidaqmx package. 
        
        The contents of the nsamples, fs, and nPxx widget assign the number of
        data points per acquisition, the sampling rate, and the number of total
        acquisitions respectively. There is a 50 ms delay implemented between 
        individual sweeps.
        
        Following each sweep, the raw data is displayed in a plot on the left
        side of the window and a mean power spectral density (across sweeps) 
        of the position is displayed in a plot on the right. The contents of 
        the sensitivity widget are used for the conversion from power in V^2
        to m^2.
        
        """
        
        ## Call variables from QLineEdit containers in GUI
        self.nsamples = int(self.ui.nsamplesinput.text())
        self.fs = int(self.ui.fsinput.text())
        self.nPxx = int(self.ui.nspectrainput.text())
        
        ## Calculate the frequency resolution of the power spectra
        self.fres = self.fs/self.nsamples
        #print('Frequency Resolution: ' + str(self.fres))
        

        ## Create a time array for plotting purposes in ms                         
        self.t = np.linspace(0, self.nsamples/self.fs,self.nsamples)
        self.t *= 1000
        
        ## Initialize arrays to contain raw data and power spectra
        self.datArr = np.zeros((self.nPxx,self.nsamples))
        self.PxxArr = np.zeros((self.nPxx,int(0.5*self.nsamples+1)))
        #self.PxxArr = np.zeros((self.nPxx,int(self.nsamples)))
        ## Acquire from nidaq device
        with nidaqmx.Task() as task:
            ## Set-up acquisition channel and configure acquisition protocol
            task.ai_channels.add_ai_voltage_chan("Dev1/ai0",
                                         terminal_config = TerminalConfiguration.BAL_DIFF,
                                         min_val = -1.0,
                                         max_val = 1.0,
                                         units = VoltageUnits.VOLTS)
            task.timing.cfg_samp_clk_timing(self.fs, samps_per_chan=self.nsamples)
            
            ## Loop during which data is acquired
            for i in range(np.size(self.datArr,axis=0)):
                
                ## Read input channel
                data = task.read(number_of_samples_per_channel=self.nsamples)
        
                ## Plot the raw data during acquisition
                self.ui.axes.cla()
                self.ui.axes.plot(self.t,data, c='b',linewidth=0.5)
                self.ui.axes.set_xlabel('Time (ms)')
                self.ui.axes.set_ylabel('Voltage (V)')
                self.ui.figure.tight_layout(pad = 0.1)
                self.ui.ainput.draw()
                
                ## Generate power spectrum
                freqs,power = signal.periodogram(data, self.fs, window='hann',
                                                 detrend ='constant',
                                                 return_onesided=True)
                self.datArr[i,:] = data
                self.PxxArr[i,:] = power
                if i > 0: 
                    self.uPxx = np.mean(self.PxxArr[:i,:],axis=0)   
                    #self.fplt = freqs[np.where(freqs >=100000)]
                    self.fplt = freqs
                    #self.uPxxplt = self.uPxx[np.where(freqs >=100000)]
                    self.uPxxplt = self.uPxx
                    self.ui.axes2.cla()
                    self.ui.axes2.plot(self.fplt,self.uPxxplt,'o',c='b',markersize=0.5)
                    self.ui.axes2.set_xlabel('Frequency (Hz)')
                    self.ui.axes2.set_ylabel('Power (V^2/Hz)')
    
                    self.ui.figure2.tight_layout(pad = 0.1)
                    self.ui.mplspectra.draw()
                sleep(0.05)
                QtCore.QCoreApplication.processEvents()
            self.freqs = freqs
        
    def onselect(self, xmin, xmax):
            indmin, indmax = np.searchsorted(self.freqs, (xmin, xmax))
            indmax = min(len(self.freqs) - 1, indmax)
            self.sub_uPxx = self.uPxx[indmin:indmax]
            self.sub_freqs = self.freqs[indmin:indmax]
        
    def data_select(self):
        self.span = SpanSelector(self.ui.axes2, self.onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.1, facecolor='grey'))
    
    def SHO(self, f,f0,a,bg,Q):
        power_app = (a*(f0**4))/((f**2-f0**2)**2 + (f*f0/Q)**2)+bg
        return power_app
    
    # function for genetic algorithm to minimize (sum of squared error)
    def sumOfSquaredError(self, parameterTuple):
        warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
        val = self.SHO(self.sub_freqs, *parameterTuple)
        return np.sum((self.sub_uPxx - val) ** 2.0)
    
    def generate_Initial_Parameters(self):
        self.parameterBounds = []
        self.parameterBounds.append([min(self.sub_freqs), max(self.sub_freqs)]) # search bounds forf0
        self.parameterBounds.append([min(self.sub_uPxx)*10**-2, max(self.sub_uPxx)]) # search bounds for a
        self.parameterBounds.append([min(self.sub_uPxx)*10**-2, max(self.sub_uPxx)]) # search bounds for b
        self.parameterBounds.append([0, 1000]) # search bounds for Q

        # "seed" the numpy random number generator for repeatable results
        self.result = differential_evolution(self.sumOfSquaredError, 
                                             self.parameterBounds, seed=3)
        return self.result.x
        
    def fit_data(self):
        # by default, differential_evolution completes by calling curve_fit() using parameter bounds
        self.geneticParameters = self.generate_Initial_Parameters()

        self.popt, self.pcov = curve_fit(self.SHO, self.sub_freqs,self.sub_uPxx, 
                                         self.geneticParameters)
        self.perr = np.sqrt(np.diag(self.pcov))
        
        #print('Fitted parameters:', self.popt)
        #print('Uncertainty:', self.perr)
        
        self.modelPredictions = self.SHO(self.freqs, *self.popt)
        
        self.ui.axes2.cla()
        self.ui.axes2.plot(self.fplt,self.uPxxplt,'o', c='b',markersize=0.5)
        self.ui.axes2.plot(self.freqs,
                           self.modelPredictions,'r-',
                           linewidth=0.5)
        self.ui.mplspectra.draw()
        
        item = self.ui.paramtable.item(1, 1)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.popt[0],2))))
        item = self.ui.paramtable.item(2, 1)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(format(self.popt[1],'.3g'))))
        item = self.ui.paramtable.item(3, 1)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.popt[3],2))))
        item = self.ui.paramtable.item(4, 1)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(format(self.popt[2],'.3g'))))
        
        item = self.ui.paramtable.item(1, 2)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.perr[0],2))))
        item = self.ui.paramtable.item(2, 2)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(format(self.perr[1],'.3g'))))
        item = self.ui.paramtable.item(3, 2)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.perr[3],2))))
        item = self.ui.paramtable.item(4, 2)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(format(self.perr[2],'.3g'))))

    def Reynolds_num(self, f0,width):
        p = 1.18
        n = 1.86e-5
    
        result = (p * f0 * width**2) / (4 * n)
        return result
    
    def hydrodynamic_circ(self,f0,width,Re):
        result = 1 + (4 * 1j * besselk(1, -1j*np.sqrt(1j*Re)))/\
        (np.sqrt(1j*Re) * besselk(0, -1j*np.sqrt(1j*Re)))
    
        return result

    def G_real(self,z):
        result = (0.91324 - 0.48274*z + 0.46842*z**2 - 0.12886*z**3\
                  + 0.044055*z**4 - 0.0035117*z**5 + 0.00069085*z**6) /\
                  (1 - 0.56964*z + 0.48690*z**2 - 0.13444*z**3\
                   + 0.045155*z**4 - 0.0035862*z**5 + 0.00069085*z**6)
        return result

    def G_imag(self,z):
        result = (-0.024134 - 0.029256*z + 0.016294*z**2 - 0.00010961*z**3\
                  + 0.000064577*z**4 - 0.00004451*z**5)/\
                  (1 - 0.59702*z + 0.55182*z**2 - 0.18357*z**3\
                   + 0.079156*z**4 - 0.014369*z**5 + 0.0028361*z**6)
        return result
    
    def hydrodynamic_corr(self,Re):
        logRe = np.log10(Re)
        result = self.G_real(logRe) + 1j*self.G_imag(logRe)
        return result

    def hydrodynamic_rect(self,f0,width):
        Re = self.Reynolds_num(f0,width)
        result = self.hydrodynamic_corr(Re) * self.hydrodynamic_circ(f0,width,Re)
    
        return result
        
    def calculate_stiffness(self):
        p = 1.18
        self.w = float(self.ui.widthedit.text()) * 1e-6
        self.l = float(self.ui.lengthedit.text()) * 1e-6
        self.dim_uncertainty = float(self.ui.uncertedit.text()) * 1e-6
        f0 = 2*np.pi*self.popt[0]
        Q = self.popt[3]

        
        self.k = 0.1906 * f0**2 * p * self.w**2 * self.l *\
                np.imag(self.hydrodynamic_rect(f0,self.w)) * Q
        item = self.ui.paramtable.item(5, 1)
        item.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.k,2))))
        self.kuncert = 2*(self.dim_uncertainty/self.w) * self.k
        uncert = self.ui.paramtable.item(5,2)
        uncert.setText(QtCore.QCoreApplication.translate("MainWindow", str(round(self.kuncert,2))))
        
    def save_raw_data(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self,'Save Raw Data',"",
                                                     "CSV Files (*.csv)")
        np.savetxt(name[0],self.datArr,delimiter=',')
        
    def save_Pv_data(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self,'Save PSD Data',"",
                                                     "CSV Files (*.csv)")
        comb = np.vstack((self.PxxArr,self.freqs))
        np.savetxt(name[0],comb,delimiter=',')
        
    
app = QtWidgets.QApplication([])
application = mywindow()
application.show()
sys.exit(app.exec_())


