# -*- coding: utf-8 -*-
"""
TEASER+
Contact:
M.Sc. Avichal Malhotra: malhotra@e3d.rwth-aachen.de
B.Sc. Maxim Shamovich: maxim.shamovich@rwth-aachen.de

www.e3d.rwth-aachen.de
Mathieustr. 30
52074 Aachen

Some GUI functions taken from CityLDT.
Please contact Avichal Malhotra or Simon Raming for any questions.
Original Repo for CityLDT: https://gitlab.e3d.rwth-aachen.de/e3d-software-tools/cityldt/cityldt.git
"""

# import of libraries
import os
import sys
import PySide2
from PySide2 import QtWidgets, QtGui
import gui_functions as gf
import tplusselection as sel
from teaser.project import Project
import simulate as sim
from teaser.logic.buildingobjects.buildingphysics.en15804lcadata import En15804LcaData


# setting environment variable for PySide2
dirname = os.path.dirname(PySide2.__file__)
plugin_path = os.path.join(dirname, 'plugins', 'platforms')
os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = plugin_path


# positions and dimensions of window
posx = 275
posy = 100
width = 650
height = 700
sizefactor = 0
sizer = True

pypath = os.path.dirname(os.path.realpath(__file__))        # path of script


class mainWindow(QtWidgets.QWidget):
    def __init__(self):
        #initiate the parent
        super(mainWindow,self).__init__()
        self.initUI()


    def initUI(self):
        global posx, posy, width, height, sizefactor, sizer

        # setup of gui / layout
        if sizer:
            posx, posy, width, height, sizefactor = gf.screenSizer(self, posx, posy, width, height, app)
            sizer = False
        gf.windowSetup(self, posx, posy, width, height, pypath, 'TEASER+')

        # Setting main layout
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)

        # Loading banner
        gf.load_banner(self, os.path.join(pypath, r'pictures\TEASER+_header.png'), sizefactor)

        #Setting Layout
        self.uGrid = QtWidgets.QGridLayout()

        self.btn_selFile = QtWidgets.QPushButton('Select file')
        self.uGrid.addWidget(self.btn_selFile, 0, 0, 1, 1)

        self.btn_selDir = QtWidgets.QPushButton('Select folder')
        self.uGrid.addWidget(self.btn_selDir, 0, 1, 1, 1)

        self.txtB_inPath = QtWidgets.QLineEdit()
        self.txtB_inPath.setPlaceholderText('Path to file or folder')
        self.txtB_inPath.setReadOnly(True)
        self.uGrid.addWidget(self.txtB_inPath, 0, 2, 1, 4)

        self.lbl_scanLoD = QtWidgets.QLabel('LoD scan progress:')
        self.uGrid.addWidget(self.lbl_scanLoD, 1, 0, 1, 1)

        self.pB_scanLoD = QtWidgets.QProgressBar(self)
        self.uGrid.addWidget(self.pB_scanLoD, 1, 1, 1, 5)

        self.vbox.addLayout(self.uGrid)

        # for selecting all or individual buildings
        self.gB_buildings = QtWidgets.QGroupBox('')
        self.vbox.addWidget(self.gB_buildings)
        # self.gB_buildings.setToolTip('When unchecked transformation will be done for all buildings in the file(s)')

        self.bGrid = QtWidgets.QGridLayout()
        self.gB_buildings.setLayout(self.bGrid)

        self.rb_allBuildings = QtWidgets.QRadioButton('Transform all buildings')
        self.bGrid.addWidget(self.rb_allBuildings, 0, 0, 1, 1)
        self.rb_allBuildings.setChecked(True)

        self.rb_selectBuildings = QtWidgets.QRadioButton('Select individual buildings')
        self.bGrid.addWidget(self.rb_selectBuildings, 0, 3, 1, 1)


        self.tbl_buildings = QtWidgets.QTableWidget()
        self.tbl_buildings.setColumnCount(4)
        self.tbl_buildings.setHorizontalHeaderLabels(['File Name', 'Name of Building', 'Level of Detail (LoD)', ''])
        self.tbl_buildings.verticalHeader().hide()
        # self.tbl_buildings.horizontalHeader().hide()
        self.tbl_buildings.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tbl_buildings.setEnabled(False)
        self.tbl_buildings.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.tbl_buildings.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.tbl_buildings.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
        self.bGrid.addWidget(self.tbl_buildings, 1, 0, 1, 6)

        # Gridbox for lower grid
        self.lGrid = QtWidgets.QGridLayout()
        self.btn_next = QtWidgets.QPushButton('TEASER Enrichment')
        self.lGrid.addWidget(self.btn_next, 0, 0, 1, 1)
        self.btn_next.setEnabled(False)

        self.btn_teasereco = QtWidgets.QPushButton('TEASEREco')
        self.lGrid.addWidget(self.btn_teasereco, 0, 2, 1, 1)
        self.btn_teasereco.setEnabled(True)

        self.btn_about = QtWidgets.QPushButton('About')
        self.lGrid.addWidget(self.btn_about, 1, 0, 1, 1)

        self.btn_reset = QtWidgets.QPushButton('Reset')
        self.lGrid.addWidget(self.btn_reset, 1, 1, 1, 1)

        self.btn_exit = QtWidgets.QPushButton('Exit')
        self.lGrid.addWidget(self.btn_exit, 1, 2, 1, 1)

        # self.btn_next = QtWidgets.QPushButton('Next')
        # self.lGrid.addWidget(self.btn_next, 0, 3, 1, 1)
        # self.btn_next.setEnabled(False)

        self.vbox.addLayout(self.lGrid)
        self.btn_selFile.clicked.connect(self.func_selectFile)
        self.btn_selDir.clicked.connect(self.func_selectDir)
        self.rb_selectBuildings.toggled.connect(self.func_selB)
        self.btn_about.clicked.connect(self.func_about)
        self.btn_reset.clicked.connect(self.func_reset)
        self.btn_exit.clicked.connect(self.func_exit)
        self.btn_next.clicked.connect(self.func_next)
        self.btn_teasereco.clicked.connect(self.func_eco)
        # self.cB_curBuilding.currentTextChanged.connect(self.func_curBuildingChanged)

        # self.gB_buildings.toggled.connect(self.func_buildingSelection)



        # setting some defaults
        self.inpPath = ''
        self.inpDir = ''
        self.expPath = ''
        self.completedLoD = 0
        # table row index to comboBox index
        self.tableDict = {}

    def func_selectFile(self):
        res = sel.select_gml(self)
        if res:
            self.inpPath = res
            self.inpDir = os.path.dirname(res)
            sel.get_files(self)
        else:
            pass

    def func_selectDir(self):
        res = sel.select_folder(self)
        if res:
            self.inpPath = res
            self.inpDir = res
            sel.get_files(self)
        else:
            pass

    def func_selB(self):
        if self.rb_selectBuildings.isChecked():
            self.tbl_buildings.setEnabled(True)
        else:
            self.tbl_buildings.setEnabled(False)

    def func_about(self):
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, about(), False)

    def func_reset(self):
        global posx, posy
        self.reset_variables()
        posx, posy = gf.dimensions(self)
        gf.next_window(self, mainWindow())

    def reset_variables(self):
        self.inpPath = ''
        self.inpDir = ''
        self.completedLoD = 0

    def func_exit(self):
        gf.close_application(self)

    def func_next(self):
        global buildingDict, selAll, inpDir
        inpDir = self.inpDir
        selAll = self.rb_allBuildings.isChecked()
        buildingDict = self.buildingDict
        gf.next_window(self, enrichment())


    def func_eco(self):
        self.btn_teasereco.setEnabled(True)
        
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, eco(), False)

    def onStateChanged(self):
        """gets called when a checkbox for a building is (un)checked to update the buildingDict"""
        ch = self.sender()
        ix = self.tbl_buildings.indexAt(ch.pos())
        self.buildingDict[ix.row()]["selected"] = ch.isChecked()
        curText = self.tbl_buildings.item(ix.row(), 1).text().split('/')[0]
        for i in range(self.tbl_buildings.rowCount()):
            if i != ix.row():
                if self.tbl_buildings.item(i, 1).text().split('/')[0] == curText:
                    self.cBoxes[i].setChecked(ch.isChecked())
                    self.buildingDict[i]["selected"] = ch.isChecked()
                    self.checkBoxChange(i, ch.isChecked())
        self.checkBoxChange(ix.row(), ch.isChecked())

    def checkBoxChange(self, row, state):
        """changes Table"""
        if state:
            colorCode = (251, 255, 0)

            # setting dummyValue to get right index of sorted list
            self.tableDict[row] = "dummyValue"
            # getting and sorting all indexes of selected buildings
            sortedList = sorted(list(self.tableDict.keys()))
            # getting the correceted Index location of the new building
            correctedIndex = sortedList.index(row) + 1
            # replacing the dummyValue
            self.tableDict[row] = correctedIndex
            # add item to comboBox
            # self.cB_curBuilding.insertItem(correctedIndex, self.buildingDict[row]["filename"] + "/" + self.buildingDict[row]["buildingname"])
            # self.btn_saveBuildingParamsAndNext.setEnabled(True)

            for y in self.tableDict:
                if y > row:
                    self.tableDict[y] = self.tableDict[y] + 1

        else:
            colorCode = (255, 255, 255)
            num = row

            # remove save building Parameters
            if row in self.buildingParamsDict:
                del self.buildingParamsDict[row]

            for y in self.tableDict:
                if y > num:
                    self.tableDict[y] = self.tableDict[y] - 1

            # remove item from comboBox
            self.cB_curBuilding.removeItem(self.tableDict[row])
            del self.tableDict[row]

            # if self.cB_curBuilding.count() == 1 or self.cB_curBuilding.currentIndex() == self.cB_curBuilding.count() - 1:
            #     self.btn_saveBuildingParamsAndNext.setEnabled(False)

        gf.setTableRowColor(self, colorCode, row)

class enrichment(QtWidgets.QWidget):
    """window for transformation options"""

    def __init__(self):
        super(enrichment, self).__init__()
        self.initUI()

    def initUI(self):
        global posx, posy, width, height, sizefactor, sizer
        if sizer:
            posx, posy, width, height, sizefactor = gf.screenSizer(self, posx, posy, width, height, app)
            sizer = False
        gf.windowSetup(self, posx, posy, width, height, pypath,
                       'CityLDT - CityGML LoD Transformation Tool - Transformation')

        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)

        gf.load_banner(self, os.path.join(pypath, r'pictures\TEASER+_header.png'), sizefactor)

        self.uGrid = QtWidgets.QGridLayout()

        self.vbox.addLayout(self.uGrid)
        ttl = 'Number of selected buildings - ' + str(len(buildingDict))
        # building parameters
        self.gB_buildingParameters = QtWidgets.QGroupBox(ttl)
        self.vbox.addWidget(self.gB_buildingParameters)
        self.vBox_forBPgB = QtWidgets.QVBoxLayout()
        self.gB_buildingParameters.setLayout(self.vBox_forBPgB)

        # building selection
        self.pGrid = QtWidgets.QGridLayout()

        self.lbl_curBuilding = QtWidgets.QLabel('Current building:')
        self.pGrid.addWidget(self.lbl_curBuilding, 0, 0, 1, 1)

        self.cB_curBuilding = QtWidgets.QComboBox()
        self.cB_curBuilding.addItems(['all (selected) buildings'])
        self.pGrid.addWidget(self.cB_curBuilding, 0, 1, 1, 2)

        presenetLoDs = []
        # adding selected buildings to the comboBox
        for key in buildingDict:
            if selAll:
                self.cB_curBuilding.insertItem(self.cB_curBuilding.count(),
                                               buildingDict[key]["filename"] + "/" + buildingDict[key]["buildingname"])
                if buildingDict[key]["values"]["LoD"] not in presenetLoDs:
                    presenetLoDs.append(buildingDict[key]["values"]["LoD"])
            elif buildingDict[key]["selected"]:
                self.cB_curBuilding.insertItem(self.cB_curBuilding.count(),
                                               buildingDict[key]["filename"] + "/" + buildingDict[key]["buildingname"])
                if buildingDict[key]["values"]["LoD"] not in presenetLoDs:
                    presenetLoDs.append(buildingDict[key]["values"]["LoD"])
            else:
                pass
        self.buildingDict = buildingDict

        # update title of groubbox according to number of buildings
        ttl = 'Number of selected buildings - ' + str(self.cB_curBuilding.count() - 1)
        self.gB_buildingParameters.setTitle(ttl)

        self.vBox_forBPgB.addLayout(self.pGrid)

        # Enrichmet properties
        self.dB_enrichment = QtWidgets.QGroupBox('Building Enrichment')
        self.vBox_forBPgB.addWidget(self.dB_enrichment)

        self.dGrid = QtWidgets.QGridLayout()
        self.dB_enrichment.setLayout(self.dGrid)

        self.lbl_year_of_construction = QtWidgets.QLabel('Year of Construction')
        self.dGrid.addWidget(self.lbl_year_of_construction, 0, 0, 1, 1)

        self.txtb_year_of_construction = QtWidgets.QLineEdit()
        # self.txtb_year_of_construction.setFixedWidth(250)
        self.txtb_year_of_construction.setPlaceholderText('Select individual building to overwrite')
        self.txtb_year_of_construction.setReadOnly(True)
        self.dGrid.addWidget(self.txtb_year_of_construction, 0, 1, 1, 1)

        self.lbl_dwelling_archetype = QtWidgets.QLabel('Select Archetype')
        self.dGrid.addWidget(self.lbl_dwelling_archetype, 0, 2, 1, 1)

        self.combobox_dwelling = QtWidgets.QComboBox()
        self.combobox_dwelling.addItems([" ", "IWU Single Family Dwelling", "TABULA Single Family House", "TABULA Multi Family House", "TABULA Terraced House", "TABULA Apartment Block"])
        self.dGrid.addWidget(self.combobox_dwelling, 0, 3, 1, 1)

        self.lbl_weather = QtWidgets.QPushButton('Select Weather file')
        self.dGrid.addWidget(self.lbl_weather, 1, 0, 1, 1)

        self.txt_weather = QtWidgets.QLineEdit('')
        self.dGrid.addWidget(self.txt_weather, 1, 1, 1, 3)

        self.lbl_enrichment = QtWidgets.QLabel('Enrich models using pre-defined archetypes')
        self.dGrid.addWidget(self.lbl_enrichment, 2, 0, 1, 2)

        self.btn_enrich = QtWidgets.QPushButton('Click to custom enrich')
        self.dGrid.addWidget(self.btn_enrich, 2, 2, 1, 2)

        self.grp_out = QtWidgets.QGroupBox('Output Options')
        self.vBox_forBPgB.addWidget(self.grp_out)

        self.mGrid = QtWidgets.QGridLayout()
        self.grp_out.setLayout(self.mGrid)

        # Adding output options
        self.checkbox_modelica = QtWidgets.QCheckBox("Modelica Model")
        self.mGrid.addWidget(self.checkbox_modelica, 0, 0, 1, 1)

        self.checkbox_withsimulation = QtWidgets.QCheckBox("With Simulation (if Dymola installed locally)")
        self.mGrid.addWidget(self.checkbox_withsimulation, 0, 1, 1, 1)

        self.checkbox_gml_ade_file = QtWidgets.QCheckBox("CityGML Energy ADE")
        self.mGrid.addWidget(self.checkbox_gml_ade_file, 0, 2, 1, 1)

        self.lbl_dymola = QtWidgets.QPushButton('Select Dymola path')
        self.mGrid.addWidget(self.lbl_dymola, 1, 0, 1, 1)
        self.lbl_dymola.setEnabled(False)

        self.txt_dymola = QtWidgets.QLineEdit('')
        self.mGrid.addWidget(self.txt_dymola, 1, 1, 1, 2)
        self.txt_dymola.setEnabled(False)

        self.grp_save = QtWidgets.QGroupBox('Save As')
        self.vBox_forBPgB.addWidget(self.grp_save)

        self.m1Grid = QtWidgets.QGridLayout()
        self.grp_save.setLayout(self.m1Grid)

        # Adding save as options
        self.checkbox_excel = QtWidgets.QCheckBox(".anything else")
        self.m1Grid.addWidget(self.checkbox_excel, 0, 0, 1, 1)

        self.checkbox_CSV = QtWidgets.QCheckBox(".csv")
        self.m1Grid.addWidget(self.checkbox_CSV, 0, 1, 1, 1)

        self.checkbox_excel = QtWidgets.QCheckBox(".xsls")
        self.m1Grid.addWidget(self.checkbox_excel, 0, 2, 1, 1)


        self.lbl_outputpath = QtWidgets.QPushButton('Select Output path')
        self.m1Grid.addWidget(self.lbl_outputpath, 1, 0, 1, 1)
        self.lbl_outputpath.setEnabled(False)

        self.txt_outputpath = QtWidgets.QLineEdit('')
        self.m1Grid.addWidget(self.txt_outputpath, 1, 1, 1, 2)
        self.txt_outputpath.setEnabled(False)

        # Gridbox for lower grid
        self.lGrid = QtWidgets.QGridLayout()

        self.btn_back = QtWidgets.QPushButton('Main Window')
        self.lGrid.addWidget(self.btn_back, 0, 0, 1, 1)

        self.btn_teasereco = QtWidgets.QPushButton('TEASEREco')
        self.lGrid.addWidget(self.btn_teasereco, 0, 1, 1, 1)

        self.btn_next = QtWidgets.QPushButton('Execute')
        self.lGrid.addWidget(self.btn_next, 0, 2, 1, 1)

        self.btn_about = QtWidgets.QPushButton('About')
        self.lGrid.addWidget(self.btn_about, 1, 0, 1, 1)

        self.btn_reset = QtWidgets.QPushButton('Reset')
        self.lGrid.addWidget(self.btn_reset, 1, 1, 1, 1)

        self.btn_exit = QtWidgets.QPushButton('Exit')
        self.lGrid.addWidget(self.btn_exit, 1, 2, 1, 1)

        self.vbox.addLayout(self.lGrid)



        # #Adding Label and button for file selection
        # self.label_selectfile = QtWidgets.QLabel("Select File")
        # self.uGrid.addWidget(self.label_selectfile, 0, 0, 1, 1)
        #
        # self.btn_selectfile = QtWidgets.QPushButton("Click to select")
        # self.uGrid.addWidget(self.btn_selectfile, 0, 1, 1, 1)
        #
        # #Adding label and combobox for building selection
        # self.label_selectbuilding = QtWidgets.QLabel("Select Building")
        # self.uGrid.addWidget(self.label_selectbuilding, 1, 0, 1, 1)
        #
        # self.dropdown_buildings = QtWidgets.QComboBox()
        # self.uGrid.addWidget(self.dropdown_buildings, 1, 1, 1, 1)
        #
        # #Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 2, 0, 1, 2)

        # #Adding Enrichment
        # self.label_enrichment = QtWidgets.QLabel("Enrichment")
        # self.uGrid.addWidget(self.label_enrichment, 3, 0, 1, 1)
        #
        # #Adding enrichment options
        # self.checkbox_gml = QtWidgets.QCheckBox("Select GML File")
        # self.uGrid.addWidget(self.checkbox_gml, 4, 0, 1, 1)
        #
        # self.checkbox_ade = QtWidgets.QCheckBox("Select CityGML + Energy ADE")
        # self.uGrid.addWidget(self.checkbox_ade, 4, 1, 1, 1)
        #
        # #Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 5, 0, 1, 2)
        #
        # #Add select dwelling
        # self.label_dwelling = QtWidgets.QLabel("Select Dwelling Archetype")
        # self.uGrid.addWidget(self.label_dwelling, 6, 0, 1, 1)
        #
        # self.combobox_dwelling = QtWidgets.QComboBox()
        # self.combobox_dwelling.addItems(["SFD", "MFD", "etc"])
        # self.uGrid.addWidget(self.combobox_dwelling, 6, 1, 1, 1)
        #
        # #Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 7, 0, 1, 2)
        #
        # #Add select weather
        # self.label_weather = QtWidgets.QLabel("Select weather file")
        # self.uGrid.addWidget(self.label_weather, 8, 0, 1, 1)
        #
        # self.btn_weather = QtWidgets.QPushButton("Click to select")
        # self.uGrid.addWidget(self.btn_weather, 8, 1, 1, 1)
        #
        # # Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 9, 0, 1, 2)
        #
        # #Add select weather
        # self.label_custom = QtWidgets.QLabel("Custom Enrichment")
        # self.uGrid.addWidget(self.label_custom, 10, 0, 1, 1)
        #
        # self.btn_custom = QtWidgets.QPushButton("Click for custom enrichment")
        # self.uGrid.addWidget(self.btn_custom, 10, 1, 1, 1)
        #
        # # Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 11, 0, 1, 2)
        #
        # # Adding output options
        # self.checkbox_withsimulation = QtWidgets.QCheckBox("With Simulation (only if Dymola installed locally)")
        # self.uGrid.addWidget(self.checkbox_withsimulation, 12, 0, 1, 1)
        #
        # self.checkbox_modelica = QtWidgets.QCheckBox("Modelica Model")
        # self.uGrid.addWidget(self.checkbox_modelica, 12, 1, 1, 1)
        #
        # self.checkbox_gml_ade_file = QtWidgets.QCheckBox("CityGML Energy ADE")
        # self.uGrid.addWidget(self.checkbox_gml_ade_file, 13, 0, 1, 1)
        #
        # self.checkbox_csv = QtWidgets.QCheckBox("CSV Results (only if Dymola installed locally)")
        # self.uGrid.addWidget(self.checkbox_ade, 13, 1, 1, 1)
        #
        # # Adding horizontal partition
        # self.uGrid.addWidget(gf.divider(), 14, 0, 1, 2)
        #
        # self.btn_execute = QtWidgets.QPushButton("Click to execute")
        # self.uGrid.addWidget(self.btn_execute, 15, 0, 1, 1)




        #Adding grid to main layout
        self.vbox.addLayout(self.uGrid)


class about(QtWidgets.QWidget):
    def __init__(self):
        super(about, self).__init__()
        self.initUI()

    def initUI(self):
        global posx, posy, width, height, sizefactor
        

        gf.windowSetup(self, posx + 10, posy + 10, width, height, pypath, 'CityBIT - About')

        # creating main layout
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)

        gf.load_banner(self, os.path.join(pypath, r'pictures\e3dHeader.png'), sizefactor)

        self.textwidget = QtWidgets.QTextBrowser()
        self.vbox.addWidget(self.textwidget)
        self.textwidget.setFontPointSize(14)
        with open(os.path.join(pypath, 'about/about.txt'), 'r') as file:
            text = file.read()
        self.textwidget.setText(text)

        self.lGrid = QtWidgets.QGridLayout()

        self.btn_repo = QtWidgets.QPushButton('Open repository')
        self.lGrid.addWidget(self.btn_repo, 0, 0, 1, 1)

        self.btn_close = QtWidgets.QPushButton('Close')
        self.lGrid.addWidget(self.btn_close, 0, 1, 1, 1)

        self.vbox.addLayout(self.lGrid)

        self.btn_repo.clicked.connect(self.open_repo)
        self.btn_close.clicked.connect(self.close_about)

    def open_repo(self):
        os.startfile('www.e3d.rwth-aachen.de')

    def close_about(self):
        self.hide()
        
class eco(QtWidgets.QWidget):
    """ Window for TEASER+eco
    """
    
    def __init__(self):          
        super(eco, self).__init__()
        self.prj = Project()
        self.prj.used_library_calc = "AixLib"
        self.initUI()
        
        self.building_groups = []
        
    def initUI(self):
        global posx, posy, width, height, sizefactor
        
        

        gf.windowSetup(self, posx + 10, posy + 10, width, height, pypath, 'Teaser+eco')
        

        # creating main layout
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)

        
        
        self.gB_buildings = QtWidgets.QGroupBox('')
        self.vbox.addWidget(self.gB_buildings)
        
        self.btn_add_bldg = QtWidgets.QPushButton('Add building')
        self.vbox.addWidget(self.btn_add_bldg)
        

        
        self.bGrid = QtWidgets.QGridLayout()
        self.gB_buildings.setLayout(self.bGrid)
        
        
        
        self.tbl_buildings = QtWidgets.QTableWidget()
        self.tbl_buildings.setColumnCount(6)
        self.tbl_buildings.setHorizontalHeaderLabels(['#', 'Name', 'type', 'Year','net leased area', 'Quantity'])
        self.tbl_buildings.verticalHeader().hide()
        # self.tbl_buildings.horizontalHeader().hide()
        
        self.tbl_buildings.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        self.tbl_buildings.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.tbl_buildings.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        

        self.bGrid.addWidget(self.tbl_buildings, 1, 0, 1, 6)
        
        """
        self.textwidget = QtWidgets.QTextBrowser()
        self.vbox.addWidget(self.textwidget)
        self.textwidget.setFontPointSize(14)
        text = "Testtext"
        self.textwidget.setText(text)
        """

        self.lGrid = QtWidgets.QGridLayout()

        self.btn_start = QtWidgets.QPushButton('Start simulation')
        self.lGrid.addWidget(self.btn_start, 0, 0, 1, 1)

        self.btn_close = QtWidgets.QPushButton('Close')
        self.lGrid.addWidget(self.btn_close, 0, 1, 1, 1)

        self.vbox.addLayout(self.lGrid)
        
        self.btn_add_bldg.clicked.connect(self.add_building)
        self.btn_start.clicked.connect(self.start)
        self.btn_close.clicked.connect(self.close)
        
    def add_building(self):
        """Function to open the addBuilding-Window
        """
        
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, addBuilding(self.prj, self), False)

    def start(self):
        """Function to open the startSimulation-Window
        """
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, startSimulation(self.prj, self), False)

    def close(self):
        """Function to close the eco-window
        """
        self.hide()

class addBuilding(QtWidgets.QWidget):
    """Window to add buildings to the TEASER+eco-project
    """
    def __init__(self, prj, parent):
        #initiate the parent
        self.prj = prj
        self.parent = parent
        super(addBuilding, self).__init__()
        self.initUI()
        

    def initUI(self):
        global posx, posy, width, height, sizefactor, sizer
        gf.windowSetup(self, posx + width+10, posy, 300, 400, pypath, 'Add building')
        
        self.only_int = QtGui.QIntValidator()
        self.only_double = QtGui.QDoubleValidator()
        
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)
        
        self.gB_parameters = QtWidgets.QGroupBox('Enrichment')
        self.vbox.addWidget(self.gB_parameters)
        
        self.pGrid = QtWidgets.QGridLayout()
        self.gB_parameters.setLayout(self.pGrid)
        
        self.lbl_method = QtWidgets.QLabel('Method:')
        self.pGrid.addWidget(self.lbl_method, 0, 1, 1, 1)
        
        self.cb_method = QtWidgets.QComboBox()
        self.pGrid.addWidget(self.cb_method, 0, 2, 1, 1)
        self.cb_method.addItems(["", 'iwu', 'urbanrenet', 'tabula_de'])
        
        self.lbl_usage = QtWidgets.QLabel('Usage:')
        self.pGrid.addWidget(self.lbl_usage, 1, 1, 1, 1)
        
        self.cb_usage = QtWidgets.QComboBox()
        self.pGrid.addWidget(self.cb_usage, 1, 2, 1, 1)
        self.cb_usage.addItems([""])
        self.cb_usage.setEnabled(False)
        
        self.lbl_year = QtWidgets.QLabel('Year of construction:')
        self.pGrid.addWidget(self.lbl_year, 2, 1, 1, 1)
        
        self.led_year = QtWidgets.QLineEdit('')
        self.pGrid.addWidget(self.led_year, 2, 2, 1, 1)
        self.led_year.setValidator(self.only_int)
        
        self.lbl_numb_flr = QtWidgets.QLabel('Number of floors:')
        self.pGrid.addWidget(self.lbl_numb_flr, 3, 1, 1, 1)
        
        self.led_numb_flr = QtWidgets.QLineEdit('')
        self.pGrid.addWidget(self.led_numb_flr, 3, 2, 1, 1)
        self.led_numb_flr.setValidator(self.only_int)
        
        self.lbl_height_flr = QtWidgets.QLabel('Height of floors:')
        self.pGrid.addWidget(self.lbl_height_flr, 4, 1, 1, 1)
        
        self.led_height_flr = QtWidgets.QLineEdit('')
        self.pGrid.addWidget(self.led_height_flr, 4, 2, 1, 1)
        self.led_height_flr.setValidator(self.only_double)
        
        self.lbl_name = QtWidgets.QLabel('Name:')
        self.pGrid.addWidget(self.lbl_name, 5, 1, 1, 1)
        
        self.led_name = QtWidgets.QLineEdit('')
        self.pGrid.addWidget(self.led_name, 5, 2, 1, 1)
        
        self.lbl_nla = QtWidgets.QLabel('Net leased area:')
        self.pGrid.addWidget(self.lbl_nla, 6, 1, 1, 1)
        
        self.led_nla = QtWidgets.QLineEdit('')
        self.pGrid.addWidget(self.led_nla, 6, 2, 1, 1)
        
        self.lbl_quantity = QtWidgets.QLabel('Quantity:')
        self.pGrid.addWidget(self.lbl_quantity, 7, 1, 1, 1)
        
        self.led_quantity = QtWidgets.QLineEdit('1')
        self.pGrid.addWidget(self.led_quantity, 7, 2, 1, 1)
        self.led_quantity.setValidator(self.only_int)
                
        self.gB_add_lca = QtWidgets.QGroupBox('Additional LCA-Data')
        self.vbox.addWidget(self.gB_add_lca)
        
        self.lBox =QtWidgets.QVBoxLayout(self)
        self.gB_add_lca.setLayout(self.lBox)
        
        self.tbl_lca = QtWidgets.QTableWidget()
        self.tbl_lca.setColumnCount(2)
        self.tbl_lca.setHorizontalHeaderLabels(['Name', 'Quantity'])
        self.tbl_lca.verticalHeader().hide()
        # self.tbl_lca.horizontalHeader().hide()
        self.tbl_lca.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tbl_lca.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.tbl_lca.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)     
        self.lBox.addWidget(self.tbl_lca)
        
        self.btn_add_lca = QtWidgets.QPushButton('Add LCA-Data')
        self.lBox.addWidget(self.btn_add_lca)
        
        self.btn_add = QtWidgets.QPushButton('Add')
        self.vbox.addWidget(self.btn_add)
        
        self.cb_method.currentTextChanged.connect(self.method_changed)
        self.btn_add.clicked.connect(self.add)
        self.btn_add_lca.clicked.connect(self.add_lca)
        
        self.additonal_lca_data = [] #list with LCA-Data IDs
        
    def method_changed(self, value):
        """Function to set suitable archetypes for the methods in the usage 
        combo box

        Parameters
        ----------
        value : str
            method (e.g. 'iwu').

        """

        self.cb_usage.clear()
        
        if value == "":
            self.cb_usage.addItem("")
            self.cb_usage.setEnabled(False)
        else:
            self.cb_usage.setEnabled(True)
            if value == "iwu":
                self.cb_usage.addItem('single_family_dwelling')
            elif value == "urbanrenet":
                self.cb_usage.addItems(['est1a', 'est1b', 'est2', 'est3', 'est4a', 'est4b', 'est5' 'est6', 'est7', 'est8a', 'est8b'])
            elif value == "tabula_de":
                self.cb_usage.addItems(["single_family_house","terraced_house","multi_family_house","apartment_block"])
                    
    def add(self):
        """Function to add the building to the project
        """
        
        method = self.cb_method.currentText()
        usage = self.cb_usage.currentText()
        year_of_construction = self.led_year.text()
        number_of_floors = self.led_numb_flr.text()
        height_of_floors = self.led_height_flr.text()
        name = self.led_name.text()
        net_leased_area = self.led_nla.text()
        building_quantity = self.led_quantity.text()
        
        if method != "" and usage != "" and year_of_construction != "" and number_of_floors != "" and height_of_floors != "" and name != "" and building_quantity != "" and net_leased_area != "":
            building_quantity = int(building_quantity)
            net_leased_area = float(net_leased_area)
        
            for i in range(building_quantity):

                if i == 0:
                    first_index = len(self.prj.buildings)              
                
                self.prj.add_residential(method = method,
                                     usage = usage,
                                     name = name,
                                     year_of_construction = year_of_construction,
                                     number_of_floors = number_of_floors,
                                     height_of_floors = height_of_floors,
                                     net_leased_area = net_leased_area)
                

                for j in range(self.tbl_lca.rowCount()):

                    self.prj.buildings[-1].add_lca_data_template(self.tbl_lca.item(j,0).text(), float(self.tbl_lca.item(j,1).text()))
                
                
                   
            if building_quantity == 1:
                index = str(first_index)
            else:
                index = f"{first_index} - {len(self.prj.buildings)-1}"
            
            rowPosition = self.parent.tbl_buildings.rowCount()
            self.parent.tbl_buildings.insertRow(rowPosition)
    
            self.parent.tbl_buildings.setItem(rowPosition , 0, QtWidgets.QTableWidgetItem(index))
            self.parent.tbl_buildings.setItem(rowPosition , 1, QtWidgets.QTableWidgetItem(name))
            self.parent.tbl_buildings.setItem(rowPosition , 2, QtWidgets.QTableWidgetItem(f"{method} {usage}"))
            self.parent.tbl_buildings.setItem(rowPosition , 3, QtWidgets.QTableWidgetItem(str(year_of_construction)))
            self.parent.tbl_buildings.setItem(rowPosition , 4, QtWidgets.QTableWidgetItem(str(net_leased_area)))
            self.parent.tbl_buildings.setItem(rowPosition , 5, QtWidgets.QTableWidgetItem(str(building_quantity)))
            
            self.parent.building_groups.append([first_index, len(self.prj.buildings)-1])
            
    def add_lca(self):
        """Function to open the addLca-Window
        """
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, addLca(self.prj, self), False)

class addLca(QtWidgets.QWidget):
    """Window to add additional LCA-data to the building.
    """
    def __init__(self,prj,parent):
        self.prj = prj
        self.parent = parent
        super(addLca, self).__init__()
        self.initUI()
        
    def initUI(self):   
        global posx, posy, width, height, sizefactor, sizer
        gf.windowSetup(self, posx + 300, posy , 300, 400, pypath, 'Add LCA-data')
        
        self.only_double = QtGui.QDoubleValidator()
        
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)
        
        self.gB_lca = QtWidgets.QGroupBox('LCA-data')
        self.vbox.addWidget(self.gB_lca)
        
        self.lGrid = QtWidgets.QGridLayout()
        self.gB_lca.setLayout(self.lGrid)
        
        self.only_double = QtGui.QDoubleValidator()
        
        self.lbl_uuid = QtWidgets.QLabel('LCA-data id')
        self.lGrid.addWidget(self.lbl_uuid, 0, 1, 1, 1)
        
        self.led_uuid = QtWidgets.QLineEdit('')
        self.lGrid.addWidget(self.led_uuid, 0, 2, 1, 1)
        
        self.lbl_amount = QtWidgets.QLabel('Amount')
        self.lGrid.addWidget(self.lbl_amount, 1, 1, 1, 1)
        
        self.led_amount = QtWidgets.QLineEdit('')
        self.lGrid.addWidget(self.led_amount, 1, 2, 1, 1)
        self.led_amount.setValidator(self.only_double)
        
        self.btn_add = QtWidgets.QPushButton('Add')
        self.vbox.addWidget(self.btn_add)
        
        self.btn_add.clicked.connect(self.add)
    
    def add(self):
        """Function to add the additional LCA-data to the building
        """
        lca_id = self.led_uuid.text()
        amount = self.led_amount.text()
        
        rowPosition = self.parent.tbl_lca.rowCount()
        self.parent.tbl_lca.insertRow(rowPosition)

        self.parent.tbl_lca.setItem(rowPosition , 0, QtWidgets.QTableWidgetItem(lca_id))
        self.parent.tbl_lca.setItem(rowPosition , 1, QtWidgets.QTableWidgetItem(amount))
        
        self.parent.tbl_lca.resizeRowsToContents()
   
    
        
    
class startSimulation(QtWidgets.QWidget):
    """Window to start the simulation
    """
    def __init__(self, prj, parent):

        self.prj = prj
        self.parent = parent
        super(startSimulation, self).__init__()
        self.initUI()
        

    def initUI(self):           
        global posx, posy, width, height, sizefactor, sizer
        gf.windowSetup(self, posx + width+10, posy+300, 500, 200, pypath, 'Start Simulation')
        
        self.dGrid = QtWidgets.QGridLayout()
        self.setLayout(self.dGrid)

        self.btn_selDir_out = QtWidgets.QPushButton('Select TEASEROurput folder')
        self.dGrid.addWidget(self.btn_selDir_out, 0, 0, 1, 1)

        self.txtB_inPath_out = QtWidgets.QLineEdit()
        self.txtB_inPath_out.setPlaceholderText('Path to file or folder')
        self.txtB_inPath_out.setReadOnly(True)
        self.dGrid.addWidget(self.txtB_inPath_out, 0, 1, 1, 4)
        
        self.btn_selDir_csv = QtWidgets.QPushButton('Select Result csv path')
        self.dGrid.addWidget(self.btn_selDir_csv, 1, 0, 1, 1)

        self.txtB_inPath_csv = QtWidgets.QLineEdit()
        self.txtB_inPath_csv.setPlaceholderText('Path to file or folder')
        self.txtB_inPath_csv.setReadOnly(True)
        self.dGrid.addWidget(self.txtB_inPath_csv, 1, 1, 1, 4)
        
        
        self.lbl_elec_lca = QtWidgets.QLabel('Select EPD for electric energy')
        self.dGrid.addWidget(self.lbl_elec_lca, 2, 0, 1, 1)
        
        self.led_elec_lca = QtWidgets.QLineEdit('')
        self.dGrid.addWidget(self.led_elec_lca, 2, 1, 1, 4)
        
        self.lbl_heat_lca = QtWidgets.QLabel('Select EPD for heating')
        self.dGrid.addWidget(self.lbl_heat_lca, 3, 0, 1, 1)
        
        self.led_heat_lca = QtWidgets.QLineEdit('')
        self.dGrid.addWidget(self.led_heat_lca, 3, 1, 1, 4)
        
        self.lbl_eff = QtWidgets.QLabel('Energy conversion efficiency')
        self.dGrid.addWidget(self.lbl_eff, 4, 0, 1, 1)
        
        self.led_eff = QtWidgets.QLineEdit('')
        self.dGrid.addWidget(self.led_eff, 4, 1, 1, 4)
        
        self.lbl_time = QtWidgets.QLabel('Time')
        self.dGrid.addWidget(self.lbl_time, 5, 0, 1, 1)
        
        self.led_time = QtWidgets.QLineEdit('50')
        self.dGrid.addWidget(self.led_time, 5, 1, 1, 4)
        
        self.btn_start = QtWidgets.QPushButton('Start')
        self.dGrid.addWidget(self.btn_start, 6, 0, 1, 1)
        
        
        
        self.btn_selDir_out.clicked.connect(self.func_selectDir_out)
        self.btn_selDir_csv.clicked.connect(self.func_selectDir_csv)
        self.btn_start.clicked.connect(self.start_simulation)
        
    def start_simulation(self):
        """Function to start the simulation
        """
        
        if self.led_time != "" and self.led_eff != "" and self.txtB_inPath_csv != "" and self.txtB_inPath_out != "" and self.led_elec_lca != "" and self.led_heat_lca != "":
        
            self.prj.calc_all_buildings()
            self.prj.export_aixlib(path = self.txtB_inPath_out.text())
            
            self.prj.period_lca_scenario = int(self.led_time.text())
            
            path1 = self.txtB_inPath_out.text().replace("/", "\\")
            path2 = self.txtB_inPath_csv.text().replace("/", "\\")
            
            sim.simulate(path = path1, prj = self.prj, loading_time = 3600, result_path = path2)
            
            lca_data_elec = En15804LcaData()
            lca_data_elec.load_lca_data_template(self.led_elec_lca.text(), self.prj.data)
            
            lca_data_heat = En15804LcaData()
            lca_data_heat.load_lca_data_template(self.led_heat_lca.text(), self.prj.data)
            
            for building in self.prj.buildings:
                
                building.calc_lca_data(False, int(self.led_time.text()))
                
                building.add_lca_data_elec(lca_data_elec)
                building.add_lca_data_heating(float(self.led_eff.text()), lca_data_heat)
                
            
        
        global posx, posy
        posx, posy = gf.dimensions(self)
        gf.next_window(self, result(self.prj, self.parent), False)
                
        
        

    
        
    def func_selectDir_out(self):
        """function to select the directory"""
        dirpath = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory")       
        self.txtB_inPath_out.setText(dirpath)
        
    
    def func_selectDir_csv(self):
        """function to select the result.csv path"""
        dirpath = QtWidgets.QFileDialog.getSaveFileName(self, 'Select .csv file')
        self.txtB_inPath_csv.setText(dirpath[0])

class result(QtWidgets.QWidget):
    """Window to display the result of the life cycle assessment"""
    def __init__(self, prj, parent):

        self.prj = prj
        self.parent = parent
        super(result, self).__init__()
        self.initUI()
        

    def initUI(self):           
        global posx, posy, width, height, sizefactor, sizer
        gf.windowSetup(self, posx + width+10, posy+300, 1500, 400, pypath, 'Result')
        
        self.vbox = QtWidgets.QVBoxLayout(self)
        self.setLayout(self.vbox)
                
        self.cb_indicator = QtWidgets.QComboBox()
        self.vbox.addWidget(self.cb_indicator)
        self.cb_indicator.addItems(["pere",
                                    "pert",
                                    "penre",
                                    "penrm",
                                    "penrt",
                                    "sm",
                                    "rsf",
                                    "nrsf",
                                    "fw",
                                    "hwd",
                                    "nhwd",
                                    "rwd",
                                    "cru",
                                    "mfr",
                                    "mer",
                                    "eee",
                                    "eet",
                                    "gwp",
                                    "odp",
                                    "pocp",
                                    "ap",
                                    "ep",
                                    "adpe",
                                    "adpf"
                                    ])
        
        self.lbl_unit = QtWidgets.QLabel('Unit:')
        self.vbox.addWidget(self.lbl_unit)
        
       
        self.tbl_lca = QtWidgets.QTableWidget()
        self.tbl_lca.setColumnCount(22)
        self.tbl_lca.setHorizontalHeaderLabels(['Building', 
                                                "amount", 
                                                "a1",
                                                "a2",
                                                "a3",
                                                "a1_a3",
                                                "a4",
                                                "a5",
                                                "b1",
                                                "b2",
                                                "b3",
                                                "b4",
                                                "b5",
                                                "b6",
                                                "b7",
                                                "c1",
                                                "c2",
                                                "c3",
                                                "c4",
                                                "d",
                                                "sum",
                                                "sum wd"
                                                ])
        self.tbl_lca.verticalHeader().hide()
        # self.tbl_lca.horizontalHeader().hide()
        self.tbl_lca.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tbl_lca.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.tbl_lca.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)     
        self.vbox.addWidget(self.tbl_lca)
              
        
        self.cGrid = QtWidgets.QGridLayout()

        self.btn_copy = QtWidgets.QPushButton('Copy csv to Clipboard')
        self.cGrid.addWidget(self.btn_copy, 0, 0, 1, 1)
        
        self.checkbox_ger = QtWidgets.QCheckBox("'German Excel'")
        self.cGrid.addWidget(self.checkbox_ger, 0, 1, 1, 1)
        
        self.vbox.addLayout(self.cGrid)
        
        
        self.select_lca("pere")
        
        self.btn_copy.clicked.connect(self.copy_to_clipboard)
        
        self.cb_indicator.currentTextChanged.connect(self.select_lca)
    
    def select_lca(self, indicator):
        
        
        while (self.tbl_lca.rowCount() > 0):
            self.tbl_lca.removeRow(0)
            
        
        for building_group in self.parent.building_groups:
            
            amount = building_group[1] - building_group[0] + 1
             
            building = self.prj.buildings[building_group[0]]
            
            lca_data = building.lca_data*amount
                
            lca_data_dict = self.lca_data_to_dict(lca_data)
            
            
            content = [building.name, amount]
            
            if indicator == "pere": 
                content.extend(lca_data_dict["pere"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "pert":
                content.extend(lca_data_dict["pert"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "penre": 
                content.extend(lca_data_dict["penre"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "penrm": 
                content.extend(lca_data_dict["penrm"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "penrt": 
                content.extend(lca_data_dict["penrt"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "sm": 
                content.extend(lca_data_dict["sm"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "rsf": 
                content.extend(lca_data_dict["rsf"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "nrsf": 
                content.extend(lca_data_dict["nrsf"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "fw": 
                content.extend(lca_data_dict["fw"])
                self.lbl_unit.setText("Unit: m^3")
            elif indicator == "hwd": 
                content.extend(lca_data_dict["hwd"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "nhwd": 
                content.extend(lca_data_dict["nhwd"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "rwd": 
                content.extend(lca_data_dict["rwd"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "cru":
                content.extend(lca_data_dict["cru"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "mfr": 
                content.extend(lca_data_dict["mfr"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "mer": 
                content.extend(lca_data_dict["mer"])
                self.lbl_unit.setText("Unit: kg")
            elif indicator == "eee": 
                content.extend(lca_data_dict["eee"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "eet": 
                content.extend(lca_data_dict["eet"])
                self.lbl_unit.setText("Unit: MJ")
            elif indicator == "gwp": 
                content.extend(lca_data_dict["gwp"])
                self.lbl_unit.setText("Unit: kg CO2 eq.")
            elif indicator == "odp": 
                content.extend(lca_data_dict["odp"])
                self.lbl_unit.setText("Unit: kg R11 eq.")
            elif indicator == "pocp": 
                content.extend(lca_data_dict["pocp"])
                self.lbl_unit.setText("Unit: kg Ethene eq.")
            elif indicator == "ap": 
                content.extend(lca_data_dict["ap"])
                self.lbl_unit.setText("Unit: kg SO2 eq.")
            elif indicator == "ep": 
                content.extend(lca_data_dict["ep"])
                self.lbl_unit.setText("Unit: kg Phosphate eq.")
            elif indicator == "adpe": 
                content.extend(lca_data_dict["adpe"])
                self.lbl_unit.setText("Unit: kg Sb eq.")
            elif indicator == "adpf": 
                content.extend(lca_data_dict["adpf"])
                self.lbl_unit.setText("Unit: MJ")
              
            self.add_lca_row(content)
                    
        
    def lca_data_to_dict(self, lca_data):
        
        result = {}
        
        result['pere'] = self.indicator_to_list(lca_data.pere)
        result['pert'] = self.indicator_to_list(lca_data.pert)
        result['penre'] = self.indicator_to_list(lca_data.penre)
        result['penrm'] = self.indicator_to_list(lca_data.penrm)
        result['penrt'] = self.indicator_to_list(lca_data.penrt)
        result['sm'] = self.indicator_to_list(lca_data.sm)
        result['rsf'] = self.indicator_to_list(lca_data.rsf)
        result['nrsf'] = self.indicator_to_list(lca_data.nrsf)
        result['fw'] = self.indicator_to_list(lca_data.fw)
        result['hwd'] = self.indicator_to_list(lca_data.hwd)
        result['nhwd'] = self.indicator_to_list(lca_data.nhwd)
        result['rwd'] = self.indicator_to_list(lca_data.rwd)
        result['cru'] = self.indicator_to_list(lca_data.cru)
        result['mfr'] = self.indicator_to_list(lca_data.mfr)
        result['mer'] = self.indicator_to_list(lca_data.mer)
        result['eee'] = self.indicator_to_list(lca_data.eee)
        result['eet'] = self.indicator_to_list(lca_data.eet)
        result['gwp'] = self.indicator_to_list(lca_data.gwp)
        result['odp'] = self.indicator_to_list(lca_data.odp)
        result['pocp'] = self.indicator_to_list(lca_data.pocp)
        result['ap'] = self.indicator_to_list(lca_data.ap)
        result['ep'] = self.indicator_to_list(lca_data.ep)
        result['adpe'] = self.indicator_to_list(lca_data.adpe)
        result['adpf'] = self.indicator_to_list(lca_data.adpf)
        
        return result
    
    def indicator_to_list(self, indicator):
        
        result = []
        
        result.append(indicator.a1)
        result.append(indicator.a2)
        result.append(indicator.a3)
        result.append(indicator.a1_a3)
        result.append(indicator.a4)
        result.append(indicator.a5)
        result.append(indicator.b1)
        result.append(indicator.b2)
        result.append(indicator.b3)
        result.append(indicator.b4)
        result.append(indicator.b5)
        result.append(indicator.b6)
        result.append(indicator.b7)
        result.append(indicator.c1)
        result.append(indicator.c2)
        result.append(indicator.c3)
        result.append(indicator.c4)
        result.append(indicator.d)
        result.append(indicator.sum_stages(False))
        result.append(indicator.sum_stages(True))

        return result

        
    
    def add_lca_row(self, row):
        
        rowPosition = self.tbl_lca.rowCount()
        self.tbl_lca.insertRow(rowPosition)

        self.tbl_lca.setItem(rowPosition , 0, QtWidgets.QTableWidgetItem(str(row[0])))
        self.tbl_lca.setItem(rowPosition , 1, QtWidgets.QTableWidgetItem(str(row[1])))
        self.tbl_lca.setItem(rowPosition , 2, QtWidgets.QTableWidgetItem(str(row[2])))
        self.tbl_lca.setItem(rowPosition , 3, QtWidgets.QTableWidgetItem(str(row[3])))
        self.tbl_lca.setItem(rowPosition , 4, QtWidgets.QTableWidgetItem(str(row[4])))
        self.tbl_lca.setItem(rowPosition , 5, QtWidgets.QTableWidgetItem(str(row[5])))
        self.tbl_lca.setItem(rowPosition , 6, QtWidgets.QTableWidgetItem(str(row[6])))
        self.tbl_lca.setItem(rowPosition , 7, QtWidgets.QTableWidgetItem(str(row[7])))
        self.tbl_lca.setItem(rowPosition , 8, QtWidgets.QTableWidgetItem(str(row[8])))
        self.tbl_lca.setItem(rowPosition , 9, QtWidgets.QTableWidgetItem(str(row[9])))
        self.tbl_lca.setItem(rowPosition , 10, QtWidgets.QTableWidgetItem(str(row[10])))
        self.tbl_lca.setItem(rowPosition , 11, QtWidgets.QTableWidgetItem(str(row[11])))
        self.tbl_lca.setItem(rowPosition , 12, QtWidgets.QTableWidgetItem(str(row[12])))
        self.tbl_lca.setItem(rowPosition , 13, QtWidgets.QTableWidgetItem(str(row[13])))
        self.tbl_lca.setItem(rowPosition , 14, QtWidgets.QTableWidgetItem(str(row[14])))
        self.tbl_lca.setItem(rowPosition , 15, QtWidgets.QTableWidgetItem(str(row[15])))
        self.tbl_lca.setItem(rowPosition , 16, QtWidgets.QTableWidgetItem(str(row[16])))
        self.tbl_lca.setItem(rowPosition , 17, QtWidgets.QTableWidgetItem(str(row[17])))
        self.tbl_lca.setItem(rowPosition , 18, QtWidgets.QTableWidgetItem(str(row[18])))
        self.tbl_lca.setItem(rowPosition , 19, QtWidgets.QTableWidgetItem(str(row[19])))
        self.tbl_lca.setItem(rowPosition , 20, QtWidgets.QTableWidgetItem(str(row[20])))
        self.tbl_lca.setItem(rowPosition , 21, QtWidgets.QTableWidgetItem(str(row[21])))
        
    def copy_to_clipboard(self):
        """Function to copy the result to clipboard
        """
        
        self.clipboard = QtGui.QClipboard()
        
        self.clipboard.setText(self.tbl_to_csv(self.checkbox_ger.isChecked()))
        
       
    def tbl_to_csv(self, german = False):
      
        
       csv = "Building,amount,a1,a2,a3,a1_a3,a4,a5,b1,b2,b3,b4,b5,b6,b7,c1,c2,c3,c4,d,sum,sum+d\n"
       
       for row in range(self.tbl_lca.rowCount()):
           for column in range(self.tbl_lca.columnCount()):
               

               csv += self.tbl_lca.item(row, column).text()
               csv += ","
                              
           csv += "\n"
           
       csv.replace("None", "")
       
       if german is True:
           csv = csv.replace(",", ";")
           csv = csv.replace(".", ",")
           
       return csv
               
               
            

        
        

if __name__ == "__main__":
    
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    
    app.setStyle('Fusion')
    widget = mainWindow()
    widget.show()
    sys.exit(app.exec_())
    
