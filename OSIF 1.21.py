# python 3 tkinter import section

from tkinter import *
import tkinter as Tkinter
import tkinter.filedialog as tkFileDialog
import tkinter.messagebox as tkMessageBox

# end python 3 tkinter import section


# python 2 Tkinter import section
'''
from Tkinter import *
import Tkinter
import tkFileDialog
import tkMessageBox
'''
# end python 2 Tkinter import section


import xlrd
import matplotlib

matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # Allows for custom positioning of sub plots in matlibplot.
import scipy.optimize

# from default python modules
import os
import re
import webbrowser
import sys

# output pretty title, version info and citation prompt.
print('\n\n\n#########################################################################')
print('python version: ' + sys.version)
print('Tkinter version: ' + str(Tkinter.TkVersion))
print('#########################################################################\n\n')
print('############################################################')
print('#####     Open Source Impedance Fitter (OSIF) v1.21    #####')
print('############################################################')
print(
    "--------------------------\nWritten by Jason Pfeilsticker for the Hydrogen fuel cell manufacturing group\nat the National Renewable Energy Lab. Feb, 2018\nV1.21 adapted Oct. 2018\n")
print(
    '\nthis program uses the matplotlib, scipy, and numpy modules and was written in python.\nIf you publish data from this program, please cite them appropriately.\n--------------------------\n\n\n')


# Main program class which is called on in line 868-ish to run the program.
class OSIF:

    def __init__(self, master):
        master.title("Open Source Impedance Fitter (OSIF) v1.21")
        master.grid()
        buttonFrame = Frame(master, pady=10, )
        InputFrame = Frame(master, padx=10)
        OutputFrame = Frame(master, padx=10)
        self.plotFrame = Frame(master, bg='blue')
        self.plotFrameToolBar = Frame(master, bg='red')

        Grid.grid_columnconfigure(buttonFrame, 0, weight=1)
        Grid.grid_rowconfigure(buttonFrame, 0, weight=1)
        Grid.grid_columnconfigure(self.plotFrame, 0, weight=1)
        Grid.grid_rowconfigure(self.plotFrame, 0, weight=1)
        Grid.grid_columnconfigure(self.plotFrameToolBar, 0, weight=1)
        Grid.grid_rowconfigure(self.plotFrameToolBar, 0, weight=1)

        buttonFrame.grid(row=0, columnspan=2)
        InputFrame.grid(row=1, column=0, sticky=N, pady=3)
        OutputFrame.grid(row=1, column=1, sticky=N, pady=3)
        self.plotFrame.grid(row=2, pady=1, padx=8, columnspan=5, sticky=N + S + E + W)
        self.plotFrameToolBar.grid(row=3, pady=1, padx=8, columnspan=5, sticky=S + W)

        self.Rmem = Param()
        self.Rcl = Param()
        self.Qdl = Param()
        self.Phi = Param()
        self.Lwire = Param()
        self.area = Param()
        self.frequencyRange = Param()
        self.loading = Param()
        self.currentDataDir = Param()
        self.currentFileName = Tkinter.StringVar(master)
        self.currentFile = NONE
        self.avgResPer = Param()
        self.activeData = Data()

        entryFont = ("Calibri", '12')
        labelFont = ("Calibri", "12")

        sdPerColumn = 5
        sdColumn = 4
        fitValueColumn = 3
        unitColumn = 2
        nonFitUnitColumn = 2
        initValueColumn = 1
        varNameColumn = 0

        Label(InputFrame, text="Initial Values", font=labelFont).grid(row=1, column=initValueColumn, sticky=W)
        Label(OutputFrame, text="Fit Values", font=labelFont).grid(row=1, column=fitValueColumn, sticky=W)
        Label(OutputFrame, text="Estimated SE", font=labelFont).grid(row=1, column=sdColumn, sticky=W)
        Label(OutputFrame, text="SE % of fit value", font=labelFont).grid(row=1, column=sdPerColumn, sticky=W)

        ################################################
        ############ INPUT INITIAL VALUES ##############
        ################################################

        Label(InputFrame, text="Rmem:", font=labelFont).grid(row=2, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[ohm*cm^2]", font=labelFont).grid(row=2, column=unitColumn, sticky=W)
        self.Rmem.IE = Entry(InputFrame, width=10, font=entryFont)
        self.Rmem.IE.grid(row=2, column=initValueColumn)
        self.Rmem.IE.insert(0, "0.03")

        Label(InputFrame, text="Rcl:", font=labelFont).grid(row=3, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[ohm*cm^2]", font=labelFont).grid(row=3, column=unitColumn, sticky=W)
        self.Rcl.IE = Entry(InputFrame, width=10, font=entryFont)
        self.Rcl.IE.grid(row=3, column=initValueColumn)
        self.Rcl.IE.insert(0, "0.1")

        Label(InputFrame, text="Qdl:", font=labelFont).grid(row=4, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[F/(cm^2*sec^phi)]", font=labelFont).grid(row=4, column=unitColumn, sticky=W)
        self.Qdl.IE = Entry(InputFrame, width=10, font=entryFont)
        self.Qdl.IE.grid(row=4, column=initValueColumn)
        self.Qdl.IE.insert(0, "2.5")

        Label(InputFrame, text="Phi:", font=labelFont).grid(row=5, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[ - ]", font=labelFont).grid(row=5, column=unitColumn, sticky=W)
        self.Phi.IE = Entry(InputFrame, width=10, font=entryFont)
        self.Phi.IE.grid(row=5, column=initValueColumn)
        self.Phi.IE.insert(0, "0.95")

        Label(InputFrame, text="Lwire:", font=labelFont).grid(row=6, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[H*cm^2]", font=labelFont).grid(row=6, column=unitColumn, sticky=W)
        self.Lwire.IE = Entry(InputFrame, width=10, font=entryFont)
        self.Lwire.IE.grid(row=6, column=initValueColumn)
        self.Lwire.IE.insert(0, "2E-5")

        Label(InputFrame, text="Cell Area:", font=labelFont).grid(row=7, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[cm^2]", font=labelFont).grid(row=7, column=nonFitUnitColumn, sticky=W)
        self.area.IE = Entry(InputFrame, width=10, font=entryFont)
        self.area.IE.grid(row=7, column=initValueColumn)
        self.area.IE.insert(0, "50")

        Label(InputFrame, text="catalyst loading:", font=labelFont).grid(row=8, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[mg/cm^2]", font=labelFont).grid(row=8, column=nonFitUnitColumn, sticky=W)
        self.loading.IE = Entry(InputFrame, width=10, font=entryFont)
        self.loading.IE.grid(row=8, column=initValueColumn)
        self.loading.IE.insert(0, "0.1")

        Label(InputFrame, text="upper Frequency bound:", font=labelFont).grid(row=9, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[Hz]", font=labelFont).grid(row=9, column=nonFitUnitColumn, sticky=W)
        self.frequencyRange.IE = Entry(InputFrame, width=10, font=entryFont)
        self.frequencyRange.IE.grid(row=9, column=initValueColumn)
        self.frequencyRange.IE.insert(0, "10000")

        Label(InputFrame, text="lower Frequency bound:", font=labelFont).grid(row=10, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[Hz]", font=labelFont).grid(row=10, column=nonFitUnitColumn, sticky=W)
        self.frequencyRange.OE = Entry(InputFrame, width=10, font=entryFont)
        self.frequencyRange.OE.grid(row=10, column=initValueColumn)
        self.frequencyRange.OE.insert(0, "1")

        ################################################
        ########### OUTPUT VALUES FROM FIT #############
        ################################################
        ioBoxWidth = 10
        self.Rmem.OE = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rmem.OE.grid(row=2, column=fitValueColumn, sticky=W)
        self.Rmem.OE.insert(0, "---")
        self.Rmem.OE.config(state='readonly')

        self.Rcl.OE = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rcl.OE.grid(row=3, column=fitValueColumn, sticky=W)
        self.Rcl.OE.insert(0, "---")
        self.Rcl.OE.config(state='readonly')

        self.Qdl.OE = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Qdl.OE.grid(row=4, column=fitValueColumn, sticky=W)
        self.Qdl.OE.insert(0, "---")
        self.Qdl.OE.config(state='readonly')

        self.Phi.OE = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Phi.OE.grid(row=5, column=fitValueColumn, sticky=W)
        self.Phi.OE.insert(0, "---")
        self.Phi.OE.config(state='readonly')

        self.Lwire.OE = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Lwire.OE.grid(row=6, column=fitValueColumn, sticky=W)
        self.Lwire.OE.insert(0, "---")
        self.Lwire.OE.config(state='readonly')

        ################################################
        ########### OUTPUT VALUE SD values #############
        ################################################

        self.Rmem.OESD = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rmem.OESD.grid(row=2, column=sdColumn, sticky=W)
        self.Rmem.OESD.insert(0, "---")
        self.Rmem.OESD.config(state='readonly')

        self.Rcl.OESD = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rcl.OESD.grid(row=3, column=sdColumn, sticky=W)
        self.Rcl.OESD.insert(0, "---")
        self.Rcl.OESD.config(state='readonly')

        self.Qdl.OESD = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Qdl.OESD.grid(row=4, column=sdColumn, sticky=W)
        self.Qdl.OESD.insert(0, "---")
        self.Qdl.OESD.config(state='readonly')

        self.Phi.OESD = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Phi.OESD.grid(row=5, column=sdColumn, sticky=W)
        self.Phi.OESD.insert(0, "---")
        self.Phi.OESD.config(state='readonly')

        self.Lwire.OESD = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Lwire.OESD.grid(row=6, column=sdColumn, sticky=W)
        self.Lwire.OESD.insert(0, "---")
        self.Lwire.OESD.config(state='readonly')

        ################################################
        ############ OUTPUT SD % OF VALUES #############
        ################################################

        self.Rmem.OESDP = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rmem.OESDP.grid(row=2, column=sdPerColumn, sticky=W)
        self.Rmem.OESDP.insert(0, "---")
        self.Rmem.OESDP.config(state='readonly')

        self.Rcl.OESDP = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Rcl.OESDP.grid(row=3, column=sdPerColumn, sticky=W)
        self.Rcl.OESDP.insert(0, "---")
        self.Rcl.OESDP.config(state='readonly')

        self.Qdl.OESDP = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Qdl.OESDP.grid(row=4, column=sdPerColumn, sticky=W)
        self.Qdl.OESDP.insert(0, "---")
        self.Qdl.OESDP.config(state='readonly')

        self.Phi.OESDP = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Phi.OESDP.grid(row=5, column=sdPerColumn, sticky=W)
        self.Phi.OESDP.insert(0, "---")
        self.Phi.OESDP.config(state='readonly')

        self.Lwire.OESDP = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.Lwire.OESDP.grid(row=6, column=sdPerColumn, sticky=W)
        self.Lwire.OESDP.insert(0, "---")
        self.Lwire.OESDP.config(state='readonly')

        ################################################
        ########## OUTPUT AVG RES |Z| VALUE ############
        ################################################

        Label(OutputFrame, text="Avg. |Z| residual % of data |Z|:", font=labelFont).grid(row=8,
                                                                                         column=(sdPerColumn - 2),
                                                                                         columnspan=2, pady=20,
                                                                                         sticky=E)
        self.avgResPer.AVGRESPER = Entry(OutputFrame, width=ioBoxWidth, font=entryFont)
        self.avgResPer.AVGRESPER.grid(row=8, column=sdPerColumn, sticky=W)
        self.avgResPer.AVGRESPER.insert(0, "---")
        self.avgResPer.AVGRESPER.config(state='readonly')

        ################################################
        #################### BUTTONS ###################
        ################################################

        self.simB = Button(buttonFrame, text="Select Data Directory", command=lambda: self.SelectDataDir())
        self.simB.grid(row=0, column=0, sticky=E)

        self.simB = Button(buttonFrame, text="Model Info.", command=self.openModelInfo)
        self.simB.grid(row=0, column=1, sticky=E)

        self.simB = Button(buttonFrame, text="Citation Info.", command=self.openCitationInfo)
        self.simB.grid(row=0, column=2, sticky=W)

        self.fitB = Button(buttonFrame, text="Fit", command=self.PerformFit)
        self.fitB.grid(row=0, column=3, sticky=E)

        self.simB = Button(buttonFrame, text="Simulate", command=self.PerformSim)
        self.simB.grid(row=0, column=4, sticky=W)

        self.simB = Button(OutputFrame, text="Save Model Data", command=self.SaveData)
        self.simB.grid(row=9, column=sdPerColumn - 2, columnspan=3, sticky=N)

        Label(buttonFrame, text="Current Directory:", font=labelFont).grid(row=1, column=0, sticky=E)
        self.currentDataDir.IE = Entry(buttonFrame, width=60, font=entryFont)
        self.currentDataDir.IE.grid(row=1, column=1, sticky=W, columnspan=2)
        self.currentDataDir.IE.insert(0, "---")
        self.currentDataDir.IE.config(state='readonly')

        Label(buttonFrame, text="Current Data File:", font=labelFont).grid(row=2, column=0, sticky=E)
        self.choices = ['Data Directory not selected']
        self.fileSelectComboBox = OptionMenu(buttonFrame, self.currentFileName, *self.choices)
        self.fileSelectComboBox.grid(row=2, column=1, sticky=EW, columnspan=2)
        self.fileSelectComboBox.config(font=entryFont)

    def openModelInfo(self):
        print(
            "--------------------------\nThe model being used in this program is Eq. 2 from the opened Fuller paper.\nThe derivation can be found in their suplimentary information. If you would\nlike to use a different model, contact Jason Pfeilsticker.\n\n")
        webbrowser.open("http://jes.ecsdl.org/content/162/6/F519.full")

    def openCitationInfo(self):
        print(
            "--------------------------\nThe Program uses the the matplotlib, scipy, and numpy modules and was written in\npython (Scientific Computing in Python). Please cite them accordingly via the opend website\n\n")
        webbrowser.open("https://www.scipy.org/citing.html")

    def SelectDataDir(self):
        newDir = tkFileDialog.askdirectory(title="Select EIS data directory") + '/'
        self.currentDataDir.IE.config(state='normal')
        self.currentDataDir.IE.delete(0, END)
        self.currentDataDir.IE.insert(0, newDir)
        self.currentDataDir.IE.config(state='readonly')

        dirList = os.listdir(newDir)
        dirList = [dataFile for dataFile in dirList if
                   ('.txt' in dataFile) | ('.xls' in dataFile) | ('.xlsx' in dataFile)]
        self.fileSelectComboBox.configure(state='normal')  # Enable drop down
        menu = self.fileSelectComboBox.children['menu']

        # Clear the menu.
        menu.delete(0, 'end')
        for file in dirList:
            # Add menu items.
            menu.add_command(label=file, command=lambda v=self.currentFileName, l=file: v.set(l))
        print('Selected Data Directory: ' + self.currentDataDir.IE.get())

    def LoadSElectedFile(self):
        # if nothing is selected, exit this method
        if (len(self.currentFileName.get()) == 0) | (self.currentFileName.get() is '---'):
            print('attempt to load on null selection')
            return

        print(
            "\n===========================================\n  loading file: " + self.currentFileName.get() + '\n===========================================\n\n')
        # clear the variables of any previous data
        self.activeData.rawFrequency = []
        self.activeData.rawzPrime = []
        self.activeData.rawZdoublePrime = []
        self.activeData.rawzMod = []
        self.activeData.rawZExperimentalComplex = []
        self.activeData.rawmodZExperimentalComplex = []

        self.activeData.dataName = self.currentFileName.get()
        self.activeData.dataNameNoExt = self.activeData.dataName[0:len(self.activeData.dataName) - 4]

        # check for different file formats to parse appropriately.
        if self.currentFileName.get().endswith('.txt'):

            self.activeData.dataNameNoExt = self.activeData.dataName[0:len(self.activeData.dataName) - 4]
            self.currentFile = open(self.currentDataDir.IE.get() + self.currentFileName.get())

            # Run through the lines in the data and parse out the variables into the above lists.
            i = 0  # number of lines in the text file (data and non data)
            dataLineString = self.currentFile.readline()
            freqCol = 0
            zPrimeCol = 1
            zDoublePrimeCol = 2
            zModCol = 3
            negateImZ = 1
            if 'Frequency' in dataLineString:
                freqCol = 1
                zPrimeCol = 2
                zDoublePrimeCol = 3
                zModCol = 4
                if ("-Z''" in dataLineString):
                    print("data loaded in has -Z'' instead of Z''; negating -Z'' column.")
                    negateImZ = -1
            k = 0
            while dataLineString:

                # For debuging file input regular expressions.

                regExTest = re.match('^\d*.\d*\t\d*.\d*\t\d*.\d*', dataLineString)
                if regExTest is not None:
                    print(
                        '----------------------------------------------------------------------------------------------------------DATA LINE BELOW')
                if regExTest is None:
                    print(
                        '----------------------------------------------------------------------------------------------------------NOT DATA LINE BELOW\n' + dataLineString)

                if (len(dataLineString) > 2) & (dataLineString[0] != '#') & (
                        re.match('^\d*.\d*\t\d*.\d*\t\d*.\d*', dataLineString) is not None) & (dataLineString != ""):
                    lineList = dataLineString.split("\t")
                    length = len(lineList)
                    last = lineList[length - 1]
                    lineList.remove(last)
                    last = last.strip()
                    lineList.insert(length - 1, last)
                    self.activeData.rawFrequency.append(float(lineList[freqCol]))
                    self.activeData.rawzPrime.append(float(lineList[zPrimeCol]))
                    self.activeData.rawZdoublePrime.append(negateImZ * float(lineList[zDoublePrimeCol]))
                    self.activeData.rawzMod.append(float(lineList[zModCol]))
                    if k == 0:
                        print('\n\tFrequency,\t\tRe(Z),\t\t\tIm(Z),\t\t\t\t|Z|')
                        k = 1
                    print(lineList[freqCol], lineList[zPrimeCol], lineList[zDoublePrimeCol], lineList[zModCol])
                dataLineString = self.currentFile.readline()
            self.currentFile.close()
            i = 0
            for real in self.activeData.rawzPrime:
                self.activeData.rawZExperimentalComplex.append((real + 1j * self.activeData.rawZdoublePrime[i]))
                self.activeData.rawmodZExperimentalComplex.append(abs(self.activeData.rawZExperimentalComplex[i]))
                i += 1

            # Change into a numpy.array type so least_squares can use it.
            self.activeData.rawFrequency = np.array(self.activeData.rawFrequency)


        # Load in default excel spread sheet output format from EIS software
        elif self.currentFileName.get().endswith('.xlsx') | self.currentFileName.get().endswith('.xls'):

            def checkForNegativeImZReturnImZ(dataRowCol):
                numRows = sheet1.col(0).__len__()
                if dataRowCol[0][3].startswith('-'):
                    i = 1
                    while i < numRows:
                        dataRowCol[i][3] = -1 * float(dataRowCol[i][3])
                        i += 1

                    newHeader = dataRowCol[0]
                    newHeader.remove(dataRowCol[0][3])
                    newHeader.insert(3, dataRowCol[0][3].strip("-"))
                    dataRowCol.remove(dataRowCol[0])
                    dataRowCol.insert(0, newHeader)
                return dataRowCol

            def sheetToListRowCol(sheet1):
                returnList = []
                i = 0
                j = 0
                numRows = sheet1.col(0).__len__()
                numCol = sheet1.row(0).__len__()
                while i < numRows:
                    tempRow = []
                    j = 0
                    while j < numCol:
                        tempRow.append(sheet1.cell(i, j).value)
                        j += 1
                    returnList.append(tempRow)
                    i += 1
                returnList = checkForNegativeImZReturnImZ(returnList)
                return returnList

            def getColDataFromData(dataRowCol, colIndex):
                i = 1
                returnCol = []
                numCol = len(dataRowCol)
                while i < numCol:
                    returnCol.append(dataRowCol[i][colIndex])
                    i += 1
                return returnCol

            xlsx = xlrd.open_workbook(self.currentDataDir.IE.get() + self.currentFileName.get())
            sheet1 = xlsx.sheet_by_index(0)
            data = sheetToListRowCol(sheet1)
            xlsx.release_resources()
            del xlsx

            self.activeData.rawFrequency = getColDataFromData(data, 1)
            self.activeData.rawzPrime = getColDataFromData(data, 2)
            self.activeData.rawZdoublePrime = getColDataFromData(data, 3)
            self.activeData.rawzMod = getColDataFromData(data, 4)
            i = 0
            print('\n\tFrequency,\t\tRe(Z),\t\t\tIm(Z),\t\t\t\t|Z|')
            for real in self.activeData.rawzPrime:
                # create things for graphing later
                self.activeData.rawZExperimentalComplex.append((real + 1j * self.activeData.rawZdoublePrime[i]))
                self.activeData.rawmodZExperimentalComplex.append(abs(self.activeData.rawZExperimentalComplex[i]))

                print(self.activeData.rawFrequency[i], self.activeData.rawzPrime[i], self.activeData.rawZdoublePrime[i],
                      self.activeData.rawzMod[i])
                i += 1

            # Change into a numpy.array type so least_squares can use it.
            self.activeData.rawFrequency = np.array(self.activeData.rawFrequency)

        print("===============================\ndone loading file\n===============================")

    def ChopFreq(self):
        tempFreq = []
        self.activeData.frequency = np.array([])
        for freq in self.activeData.rawFrequency:
            if (freq > float(self.frequencyRange.OE.get())) & (freq < float(self.frequencyRange.IE.get())):
                tempFreq.append(freq)

        # chop the data to the frequency range specified in set up
        self.activeData.frequency = np.array(tempFreq)
        minIndex = self.activeData.rawFrequency.tolist().index(self.activeData.frequency[0])
        maxIndex = self.activeData.rawFrequency.tolist().index(
            self.activeData.frequency[self.activeData.frequency.shape[0] - 1])

        self.activeData.zPrime = self.activeData.rawzPrime[minIndex:maxIndex + 1]
        self.activeData.ZdoublePrime = self.activeData.rawZdoublePrime[minIndex:maxIndex + 1]
        self.activeData.zMod = self.activeData.rawzMod[minIndex:maxIndex + 1]
        self.activeData.modZExperimentalComplex = self.activeData.rawmodZExperimentalComplex[minIndex:maxIndex + 1]

    def PerformSim(self):
        self.LoadSElectedFile()
        if len(self.activeData.rawzPrime) == 0:
            tkMessageBox.showinfo("Error!", "No data file loaded\nor data is in incorrect format")
            return

        else:
            self.ChopFreq()
            ###### /float(self.area.IE.get())      Rmem
            params = [float(self.Lwire.IE.get()) / float(self.area.IE.get()),
                      float(self.Rmem.IE.get()) / float(self.area.IE.get()),
                      float(self.Rcl.IE.get()) / float(self.area.IE.get()),
                      float(self.Qdl.IE.get()),
                      float(self.Phi.IE.get())]

            self.CreateFigures(params, 'sim')

            simResiduals = self.funcCost(params)
            self.resPercentData = np.sum(simResiduals / self.activeData.zMod * 100) / len(simResiduals)
            self.avgResPer.AVGRESPER.config(state='normal')
            self.avgResPer.AVGRESPER.delete(0, END)
            self.avgResPer.AVGRESPER.insert(0, '%5.4f' % self.resPercentData)
            self.avgResPer.AVGRESPER.config(state='readonly')


    def PerformFit(self):
        self.LoadSElectedFile()
        print('\n\n\n\n' + 'Sample: ' + self.currentFileName.get() + '\n')
        if len(self.activeData.rawzPrime) == 0:
            tkMessageBox.showinfo("Error!", "No data file loaded\nor data is in incorrect format")
            return

        else:
            ### /float(self.area.IE.get())   Rmem
            self.ChopFreq()
            params = [float(self.Lwire.IE.get()) / float(self.area.IE.get()),
                      float(self.Rmem.IE.get()) / float(self.area.IE.get()),
                      float(self.Rcl.IE.get()) / float(self.area.IE.get()),
                      float(self.Qdl.IE.get()),
                      float(self.Phi.IE.get())]

            # Perform the fitting using least_squares with the TRF method and max function calls of 10000.
            finalOutput = scipy.optimize.least_squares(self.funcCost, params, max_nfev=50000, method='trf', xtol=1e-11,
                                                       ftol=1e-11, gtol=1e-11, verbose=1)
            self.finalParams = finalOutput.x

            # Estimate variance of parameters based on Gauss-Newton approximation of the Hessian of the cost function. See: (https://www8.cs.umu.se/kurser/5DA001/HT07/lectures/lsq-handouts.pdf)
            # basically Covariance matrix = inverse(Jacob^T*Jacob)*meanSquaredError, where Jacob^T*Jacob is the first order estimate for the hessian. The square root of the diagonal elements (c_ii) of Cov are the variances of the parameter b_i

            # sigma squared estimate also called s^2 sometimes = chi^2 reduced
            sigmaSquared = (np.matmul(np.transpose(finalOutput.fun), finalOutput.fun)) / (
                    finalOutput.fun.shape[0] - self.finalParams.size)


            # Plot residuals for error checking. Note: cancels update of the normal plots in the normal UI.
            # plt.figure(100)
            # plt.plot(finalOutput.fun)
            # plt.show()

            Jacob = finalOutput.jac
            estVars = np.matrix.diagonal(
                np.linalg.inv(Jacob.T.dot(Jacob)) * (np.matmul(np.transpose(finalOutput.fun), finalOutput.fun)) / (
                            finalOutput.fun.shape[0] - self.finalParams.size))
            # estVars = sigmaSquared*np.matrix.diagonal(np.linalg.inv(np.matmul(np.matrix.transpose(finalOutput.jac),finalOutput.jac)))

            # estimated errors in parameters
            self.standardDeviation = np.sqrt(estVars)

            # taking chi^2  = sum((residuals^T)*(residuals)) = sum (r_i)^2
            # print('Sum res^2 = ' + str(np.sum((np.matrix.transpose(finalOutput.fun).__mul__(finalOutput.fun)))))

            self.L2NormOfRes = np.sqrt(np.sum(pow(finalOutput.fun, 2)))

            # print('\nNormalized grad = grad/|grad| = ' + str(finalOutput.grad / np.sqrt(np.sum(pow(finalOutput.grad, 2)))))

            self.resPercentData = np.sum(finalOutput.fun / self.activeData.zMod * 100) / finalOutput.fun.shape[0]

            print('\nFit to: ' + self.activeData.dataNameNoExt)

            self.percentSigma = self.standardDeviation / self.finalParams * 100

            self.fitOutPutString = '#\n#\n#\t\t\t\t\t\t\t   Fit values\t\t\t~std Error\t\t\t ~std Error %% of value\n#\n#\tRmem  [ohm*cm^2] \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tRcl   [ohm*cm^2] \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f' \
                                   '\n#\tQdl   [F/(cm^2*sec^phi)]  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tphi   [ ]  \t\t\t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tLwire [H*cm^2] \t\t\t  = %.4e\t\t\t%.3e\t\t\t\t%8.2f' \
                                   '\n#\n#\tQdl/mgpt = %5.6f\n#\tL2 norm of res = %10.8f [ohm*cm^2]' % \
                                   (float(self.finalParams[1]) * float(self.area.IE.get()),float(self.standardDeviation[1]) * float(self.area.IE.get()), self.percentSigma[1],
                                    float(self.finalParams[2]) * float(self.area.IE.get()),float(self.standardDeviation[2]) * float(self.area.IE.get()), self.percentSigma[2],
                                    float(self.finalParams[3]), float(self.standardDeviation[3]), self.percentSigma[3],
                                    float(self.finalParams[4]), float(self.standardDeviation[4]), self.percentSigma[4],
                                    float(self.finalParams[0]) * float(self.area.IE.get()),float(self.standardDeviation[0]) * float(self.area.IE.get()), self.percentSigma[0],
                                    float(self.finalParams[3]) / (float(self.area.IE.get()) * float(self.loading.IE.get())),
                                    float(self.L2NormOfRes)*float(self.area.IE.get()))

            print(self.fitOutPutString)
            self.realFinalModel = self.funcreal(self.finalParams)
            self.imagFinalModel = self.funcImg(self.finalParams)
            self.zModFinalModel = self.funcAbs(self.finalParams)

            self.AvgRealResPer = np.sum(abs(
                np.array(abs(self.realFinalModel - self.activeData.zPrime)) / np.array(self.activeData.zPrime)) * 100) / \
                                 self.realFinalModel.shape[0]
            self.AvgImagResPer = np.sum(abs(
                np.array(abs(self.imagFinalModel - self.activeData.ZdoublePrime)) / np.array(
                    self.activeData.ZdoublePrime)) * 100) / \
                                 self.imagFinalModel.shape[0]
            # print('\nAvgRealResPer = ' + str(self.AvgRealResPer) + '\nAvgImagResPer = ' + str(self.AvgImagResPer))

            self.Lwire.OE.config(state='normal')
            self.Lwire.OE.delete(0, END)
            self.Lwire.OE.insert(0, '%5.8f' % (float(self.finalParams[0])*float(self.area.IE.get())))
            self.Lwire.OE.config(state='readonly')

            self.Rmem.OE.config(state='normal')
            self.Rmem.OE.delete(0, END)
            self.Rmem.OE.insert(0, '%5.8f' % (float(self.finalParams[1]) * float(self.area.IE.get())))
            self.Rmem.OE.config(state='readonly')

            self.Rcl.OE.config(state='normal')
            self.Rcl.OE.delete(0, END)
            self.Rcl.OE.insert(0, '%5.8f' % (self.finalParams[2] * float(self.area.IE.get())))
            self.Rcl.OE.config(state='readonly')

            self.Qdl.OE.config(state='normal')
            self.Qdl.OE.delete(0, END)
            self.Qdl.OE.insert(0, '%5.8f' % self.finalParams[3])
            self.Qdl.OE.config(state='readonly')

            self.Phi.OE.config(state='normal')
            self.Phi.OE.delete(0, END)
            self.Phi.OE.insert(0, '%5.8f' % self.finalParams[4])
            self.Phi.OE.config(state='readonly')

            self.Lwire.OESD.config(state='normal')
            self.Lwire.OESD.delete(0, END)
            self.Lwire.OESD.insert(0, '%5.8f' % (self.standardDeviation[0] * float(self.area.IE.get())))
            self.Lwire.OESD.config(state='readonly')

            self.Rmem.OESD.config(state='normal')
            self.Rmem.OESD.delete(0, END)
            self.Rmem.OESD.insert(0, '%5.8f' % ((float(self.standardDeviation[1]) * float(self.area.IE.get()))))
            self.Rmem.OESD.config(state='readonly')

            self.Rcl.OESD.config(state='normal')
            self.Rcl.OESD.delete(0, END)
            self.Rcl.OESD.insert(0, '%5.8f' % (self.standardDeviation[2] * float(self.area.IE.get())))
            self.Rcl.OESD.config(state='readonly')

            self.Qdl.OESD.config(state='normal')
            self.Qdl.OESD.delete(0, END)
            self.Qdl.OESD.insert(0, '%5.8f' % self.standardDeviation[3])
            self.Qdl.OESD.config(state='readonly')

            self.Phi.OESD.config(state='normal')
            self.Phi.OESD.delete(0, END)
            self.Phi.OESD.insert(0, '%5.8f' % self.standardDeviation[4])
            self.Phi.OESD.config(state='readonly')

            self.Lwire.OESDP.config(state='normal')
            self.Lwire.OESDP.delete(0, END)
            self.Lwire.OESDP.insert(0, '%5.4f' % (self.percentSigma[0]))
            self.Lwire.OESDP.config(state='readonly')

            self.Rmem.OESDP.config(state='normal')
            self.Rmem.OESDP.delete(0, END)
            self.Rmem.OESDP.insert(0, '%5.4f' % (self.percentSigma[1]))
            self.Rmem.OESDP.config(state='readonly')

            self.Rcl.OESDP.config(state='normal')
            self.Rcl.OESDP.delete(0, END)
            self.Rcl.OESDP.insert(0, '%5.4f' % (self.percentSigma[2]))
            self.Rcl.OESDP.config(state='readonly')

            self.Qdl.OESDP.config(state='normal')
            self.Qdl.OESDP.delete(0, END)
            self.Qdl.OESDP.insert(0, '%5.4f' % self.percentSigma[3])
            self.Qdl.OESDP.config(state='readonly')

            self.Phi.OESDP.config(state='normal')
            self.Phi.OESDP.delete(0, END)
            self.Phi.OESDP.insert(0, '%5.4f' % self.percentSigma[4])
            self.Phi.OESDP.config(state='readonly')

            self.avgResPer.AVGRESPER.config(state='normal')
            self.avgResPer.AVGRESPER.delete(0, END)
            self.avgResPer.AVGRESPER.insert(0, '%5.4f' % self.resPercentData)
            self.avgResPer.AVGRESPER.config(state='readonly')

            self.CreateFigures(self.finalParams, 'fit')

    def CreateFigures(self, params, fitOrSim):
        if fitOrSim == 'fit':
            graphLabel = 'Full complex fit: '
        elif fitOrSim == 'sim':
            graphLabel = 'Simulated using: '
        else:
            graphLabel = ''

        plt.close('all')
        # make layout for graphs
        gs0 = gridspec.GridSpec(1, 2)
        gs00 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[0])
        gs01 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[1])
        f = plt.figure(1, figsize=[8, 3.5], tight_layout='true')

        ###########################################################
        ####          PLOT NYQUIST COMBINED FITTINGS           ####
        ###########################################################

        nyGraph = plt.Subplot(f, gs01[:, :])
        f.add_subplot(nyGraph)
        nyGraph.plot(self.activeData.zPrime, self.activeData.ZdoublePrime, 'bo', ls='--', markersize=2, linewidth=1,
                     label='data: ' + self.activeData.dataNameNoExt)
        nyGraph.plot(self.funcreal(params), self.funcImg(params), 'ro', markersize=2,
                     label='\n%s\nLwire=%.5e\nRmem=%5.8f\nRcl=%5.8f\nQdl=%5.5f\nphi=%5.5f' % (
                     graphLabel, params[0]*float(self.area.IE.get()), params[1]*float(self.area.IE.get()), params[2]*float(self.area.IE.get()), params[3], params[4]))
        plt.gca().invert_yaxis()
        plt.xticks(rotation=20)

        plt.xlabel('Re(Z)')
        plt.ylabel('Im(Z)')
        plt.legend(loc=2, fontsize=6)

        ###########################################################
        ####         PLOT |Z| vs w COMBINED FITTINGS           ####
        ###########################################################

        modZgraph = plt.Subplot(f, gs00[-3, :3])
        f.add_subplot(modZgraph)
        modZgraph.plot(self.activeData.frequency, self.activeData.modZExperimentalComplex, 'bo', ls='--', markersize=2,
                       linewidth=1)
        modZgraph.plot(self.activeData.frequency, self.funcAbs(params), 'ro', markersize=2)
        plt.ylabel('|Z|')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')

        ###########################################################
        ####                   PLOT Im(Z)                      ####
        ###########################################################

        imZgraph = plt.Subplot(f, gs00[-2, :3])
        f.add_subplot(imZgraph)
        imZgraph.plot(self.activeData.frequency, self.activeData.ZdoublePrime, 'bo', ls='--', markersize=2, linewidth=1)
        imZgraph.plot(self.activeData.frequency, (self.funcImg(params)), 'ro', markersize=2)
        plt.ylabel('Im(Z)')
        plt.gca().set_yscale('linear')
        plt.gca().set_xscale('log')
        plt.gca().set_xticks([])

        ###########################################################
        ####                   PLOT Re(Z)                      ####
        ###########################################################

        reZgraph = plt.Subplot(f, gs00[-1, :3])
        f.add_subplot(reZgraph)
        reZgraph.plot(self.activeData.frequency, self.activeData.zPrime, 'bo', ls='--', markersize=2, linewidth=1)
        reZgraph.plot(self.activeData.frequency, self.funcreal(params), 'ro', markersize=2)
        plt.xlabel('frequency')
        plt.ylabel('Re(Z)')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')

        ###########################################################
        ####              Draw figure in Tkinter               ####
        ###########################################################
        for widget in self.plotFrame.winfo_children():
            widget.destroy()
        for widget in self.plotFrameToolBar.winfo_children():
            widget.destroy()

        dataPlot = FigureCanvasTkAgg(f, master=self.plotFrame)
        dataPlot.close_event()
        dataPlot.draw()
        dataPlot.get_tk_widget().grid(row=0, sticky=N + S + E + W, )
        toolbar = NavigationToolbar2TkAgg(dataPlot, self.plotFrameToolBar)
        toolbar.update()
        dataPlot._tkcanvas.grid(row=0, sticky=W + S)

        print('done with plotting')


    def KILLALL(self):
        for widget in self.plotFrame.winfo_children():
            widget.destroy()
        for widget in self.plotFrameToolBar.winfo_children():
            widget.destroy()
        print("\n\nAll plots killed. \nHave a nice day!")

    def SaveData(self):
        if (len(self.currentFileName.get()) == 0) | (self.currentFileName.get() is '---'):
            print('no data loaded')
            tkMessageBox.showinfo("Error!", "No data file loaded")
            # nothing selected to load
            return

        dataOutFile = open(self.currentDataDir.IE.get() + self.activeData.dataNameNoExt + '_fit.txt', "w+")
        i = 0
        dataOutFile.write('#Fitted model at fitting frequencies:\n#Frequency\t\tRe(Z)\t\t\tIm(Z)\t\t\t|Z|\n')
        for real in self.realFinalModel:
            dataOutFile.write(
                str(self.activeData.frequency[i]) + '\t' + str(self.realFinalModel[i]) + '\t' + str(
                    self.imagFinalModel[i]) + '\t' + str(
                    self.zModFinalModel[i]) + '\n')
            i += 1

        dataOutFile.write(
            '#\n#\n#\t\t\t\t   Fit values\t\t\t~std dev\t\t\t   ~stdDev %% of value\n#\n#\tRmem  [ohm*cm^2] \t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tRcl   [ohm*cm^2] \t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f' \
            '\n#\tQdl   [F/(cm^2*sec^phi)]  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tphi   [ ]  \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tLwire [H*cm^2] \t\t  = %.4e\t\t\t%.3e\t\t\t\t%8.2f' \
            '\n#\n#\tQdl/mgpt = %5.6f\n#\tL2 norm of res = %10.8f' % (
                float(self.finalParams[1]) * float(self.area.IE.get()),
                float(self.standardDeviation[1]) * float(self.area.IE.get()), self.percentSigma[1],
                float(self.finalParams[2]) * float(self.area.IE.get()),
                float(self.standardDeviation[2]) * float(self.area.IE.get()), self.percentSigma[2],
                float(self.finalParams[3]), float(self.standardDeviation[3]), self.percentSigma[3],
                float(self.finalParams[4]), float(self.standardDeviation[4]), self.percentSigma[4],
                float(self.finalParams[0]) * float(self.area.IE.get()),
                float(self.standardDeviation[0]) * float(self.area.IE.get()), self.percentSigma[0],
                float(self.finalParams[3]) / (float(self.area.IE.get()) * float(self.loading.IE.get())),
                float(self.L2NormOfRes)))
        dataOutFile.write('\n#\tAvg. |Z| residual % WRT to data |Z| = ' + str(self.resPercentData))

        dataOutFile.close()
        print("saved data in: " + self.currentDataDir.IE.get() + self.activeData.dataNameNoExt + '_fit.txt')

    # Because built in coth(x) cant deal with complex numbers because exp(x) cant deal with them, but pow(x,y) can.
    def JPcoth(self, x):
        return (pow(np.e, x) + pow(np.e, -x)) / (pow(np.e, x) - pow(np.e, -x))

    # Minimizing this function results in fitting the real and complex parts of the impedance at the same time.
    def funcCost(self, params):
        return np.array(np.sqrt(pow((self.funcreal(params) - self.activeData.zPrime), 2) + pow(
            (self.funcImg(params) - self.activeData.ZdoublePrime), 2)))

    # Define the functions of the model (equation 2 from "A Physics-Based Impedance Model of Proton Exchange Membrane Fuel
    # Cells Exhibiting Low-Frequency Inductive Loops")
    def funcAbs(self, param):
        return abs(self.funcreal(param) + 1j * self.funcImg(param))

    def funcreal(self, param):
        return np.real(1j * param[0] * self.activeData.frequency + param[1] + pow(
            (param[2] / (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))), 0.5) * self.JPcoth(
            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))

    def funcImg(self, param):
        return np.imag(1j * param[0] * self.activeData.frequency + param[1] + pow(
            (param[2] / (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))), 0.5) * self.JPcoth(
            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))


class Param():

    def __init__(self):
        self.IE = Entry()
        self.OE = Entry()
        self.OESD = Entry()
        self.OESDP = Entry()
        self.AVGRESPER = Entry()


class Data():

    def __init__(self):
        self.dataName = ''
        self.dataNameNoExt = ''

        self.zPrime = []
        self.ZdoublePrime = []
        self.zMod = []
        self.modZExperimentalComplex = []
        self.frequency = np.array([])

        self.rawzPrime = []
        self.rawZdoublePrime = []
        self.rawzMod = []
        self.rawmodZExperimentalComplex = []
        self.rawFrequency = []


def on_closing():
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
        app.KILLALL()
        root.destroy()
        os._exit(0)


root = Tk()
app = OSIF(root)
root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()