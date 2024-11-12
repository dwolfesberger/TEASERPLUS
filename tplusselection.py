# import of libraries
from PySide2 import QtWidgets, QtGui
import os
import lxml.etree as ET
import glob
import math

# import of functions
import gui_functions as gf

def select_gml(self):
    """func to select file"""
    tup = QtWidgets.QFileDialog.getOpenFileName(self, 'Select .gml or .xml file', self.tr("*.gml;*.xml"))
    path = tup[0]
    if path.endswith('.gml') or path.endswith('.xml'):
        self.txtB_inPath.setText(path)
        return path
    else:
        self.txtB_inPath.setText('')
        gf.messageBox(self, "Important", "Valid File not selected")
        return 0



def select_folder(self):
    """func to select folder"""
    dirpath = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory")
    if dirpath:
        self.txtB_inPath.setText(dirpath)
        return dirpath
    else:
        self.txtB_inPath.setText('')
        gf.messageBox(self, "Important", "Valid Folder not selected")
        return 0


def get_files(self):
    """func to loop through all the files and buildings and add them to the table widget for selection"""
    #
    # function to reset the table IMPORTANT
    #
    resultsDict = {}
    if os.path.isfile(self.inpPath):
        # case for single file
        resultsDict[os.path.basename(self.inpPath)] = get_lods(self.inpPath)
        gf.progressLoD(self, 100)
        pass
    elif os.path.isdir(self.inpPath):
        # case for multiple files
        filenames = glob.glob(os.path.join(self.inpPath, "*.gml")) + glob.glob(os.path.join(self.inpPath, "*.xml"))
        for i, filename in enumerate(filenames):
            resultsDict[os.path.basename(filename)] = get_lods(filename)
            gf.progressLoD(self, (i + 1) / len(filenames) * 100)
        pass
    else:
        gf.messageBox(self, "ERROR!", "Input path is neither file or directory.\nPlease reselect input data.")

    display_file_lod(self, resultsDict)
    self.btn_next.setEnabled(True)


def get_lods(filename):
    """gets all files in a building"""
    # parsing file
    print("parsing", filename)
    tree = ET.parse(filename)
    root = tree.getroot()
    nss = root.nsmap

    buildings = {}

    # check for building parts
    # one analysis function for all building parameters:
    # bHeight, rHeight, rHeading, rType, bFunction, YOC, SAG, SBG
    # look at CITYBIT you dumb idiot

    # getting all buildings in file
    buildings_in_file = root.findall('cityObjectMember/bldg:Building', nss)

    # iterating all buildings
    for building_E in buildings_in_file:
        buildingName = building_E.attrib['{http://www.opengis.net/gml}id']
        info = get_info_from_building(building_E, nss)
        if info != {}:
            buildings[buildingName] = info
        else:
            # no ground coordinates or LoD found -> can't work with building
            pass
        bps_in_bldg = building_E.findall('./bldg:consistsOfBuildingPart', nss)
        for co_bp_E in bps_in_bldg:
            bp_E = co_bp_E.find('bldg:BuildingPart', nss)
            buildingParIDJoinded = buildingName + '/' + bp_E.attrib['{http://www.opengis.net/gml}id']
            info = get_info_from_building(bp_E, nss)
            if info != {}:
                if info["bFunction"] == 'N/D':
                    buildingFunction_E = building_E.find('bldg:function', nss)
                    if buildingFunction_E != None:
                        info["bFunction"] = buildingFunction_E.text
                buildings[buildingParIDJoinded] = info
            else:
                # no ground coordinates or LoD found -> can't work with building
                pass

    return buildings


def display_file_lod(self, filesDict):
    """adds results from get_files to table"""
    self.tbl_buildings.setRowCount(0)
    self.tbl_buildings.horizontalHeader().show()
    self.cBoxes = []

    self.buildingDict = {}

    # iterating over files
    for filename in filesDict:
        buildings = filesDict[filename]

        # checking if buildings have been found in file
        if buildings == {}:
            continue

        for i, entry in enumerate(buildings):
            rowCount = self.tbl_buildings.rowCount()
            self.tbl_buildings.insertRow(rowCount)
            if i == 0:
                # with filename
                newItem = QtWidgets.QTableWidgetItem(str(filename))
                self.tbl_buildings.setItem(rowCount, 0, newItem)
            else:
                # without filename
                newItem = QtWidgets.QTableWidgetItem("")
                self.tbl_buildings.setItem(rowCount, 0, newItem)
                pass

            newItem = QtWidgets.QTableWidgetItem(str(entry))
            self.tbl_buildings.setItem(rowCount, 1, newItem)

            newItem = QtWidgets.QTableWidgetItem(str(buildings[entry]["LoD"]))
            self.tbl_buildings.setItem(rowCount, 2, newItem)

            self.cBoxes.append(QtWidgets.QCheckBox(parent=self.tbl_buildings))
            self.cBoxes[-1].clicked.connect(self.onStateChanged)
            self.tbl_buildings.setCellWidget(rowCount, 3, self.cBoxes[-1])

            self.buildingDict[rowCount] = {"filename": filename, 'buildingname': entry, 'values': buildings[entry],
                                           "selected": False}

    self.tbl_buildings.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
    self.tbl_buildings.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
    self.tbl_buildings.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
    return


def get_info_from_building(element, nss):
    """gathers necessary info on building"""
    # bHeight, rHeight, rHeading, rType, bFunction, YOC, SAG, SBG
    data = {}
    gS_list = getGroundSurfaceCoorOfBuild(element, nss)
    # getting coordinates of groundSurface of the building
    if gS_list == '':
        # no geometry found -> skipping building
        return {}
    else:
        # found geometry of building -> can continue
        pass

    lod = get_lod(element, nss)
    if lod == -1:
        # lod is not defined -> can't continue with building
        return {}
    else:
        # found LoD -> can continue
        data["LoD"] = lod

    # getting roof Type
    rootType_E = element.find('bldg:roofType', nss)
    if rootType_E != None:
        data["rType"] = rootType_E.text
    else:
        data["rType"] = 'N/D'

    # getting building function
    buildingFunction_E = element.find('bldg:function', nss)
    if buildingFunction_E != None:
        data["bFunction"] = buildingFunction_E.text
    else:
        data["bFunction"] = 'N/D'


    # getting year of construction
    yoc_E = element.find('bldg:yearOfConstruction', nss)
    if yoc_E != None:
        data["YOC"] = yoc_E.text
    else:
        data["YOC"] = 'N/D'

    return data


def getGroundSurfaceCoorOfBuild(element, nss):
    """returns the ground surface coor form element"""

    # LoD0
    for tagName in ['bldg:lod0FootPrint', 'bldg:lod0RoofEdge']:
        LoD_zero_E = element.find(tagName, nss)
        if LoD_zero_E != None:
            posList_E = LoD_zero_E.find('.//gml:posList', nss)

            if posList_E != None:
                return get_3dPosList_from_str(posList_E.text)

            else:  # case hamburg lod2 2020
                pos_Es = LoD_zero_E.findall('.//gml:pos', nss)
                polygon = []
                for pos_E in pos_Es:
                    polygon.append(pos_E.text)
                polyStr = ' '.join(polygon)
                return get_3dPosList_from_str(polyStr)

    groundSurface_E = element.find('bldg:boundedBy/bldg:GroundSurface', nss)
    if groundSurface_E != None:
        posList_E = groundSurface_E.find('.//gml:posList', nss)  # searching for list of coordinates

        if posList_E != None:  # case aachen lod2
            return get_3dPosList_from_str(posList_E.text)

        else:  # case hamburg lod2 2020
            pos_Es = groundSurface_E.findall('.//gml:pos', nss)
            polygon = []
            for pos_E in pos_Es:
                polygon.append(pos_E.text)
            polyStr = ' '.join(polygon)
            return get_3dPosList_from_str(polyStr)

    #  checking if no groundSurface element has been found
    else:  # case for lod1 files
        geometry = element.find('bldg:lod1Solid', nss)
        if geometry != None:
            poly_Es = geometry.findall('.//gml:Polygon', nss)
            all_poylgons = []
            for poly_E in poly_Es:
                polygon = []
                posList_E = element.find('.//gml:posList', nss)  # searching for list of coordinates
                if posList_E != None:
                    polyStr = posList_E.text
                else:
                    pos_Es = poly_E.findall('.//gml:pos', nss)  # searching for individual coordinates in polygon
                    for pos_E in pos_Es:
                        polygon.append(pos_E.text)
                    polyStr = ' '.join(polygon)
                coor_list = get_3dPosList_from_str(polyStr)
                all_poylgons.append(coor_list)

            # to get the groundSurface polygon, the average height of each polygon is calculated and the polygon with the lowest average height is considered the groundsurface
            averages = []
            for polygon in all_poylgons:
                # need to get polygon with lowest z coordinate here
                average = 0
                for i in range(len(polygon) - 1):
                    average -= - polygon[i][2]
                averages.append(average / (len(polygon) - 1))

            return all_poylgons[averages.index(min(averages))]
        else:
            return ''


def get_lod(element, nss):
    """returns the first LoD found in an building or buildingPart"""
    lodFlags = {'bldg:lod0FootPrint': 0, 'bldg:lod1Solid': 1, 'bldg:lod2Solid': 2, 'bldg:lod3MultiSurface': 3,
                'bldg:lod4MultiSurface': 4}
    for flag in lodFlags:
        if element.find('./' + flag, nss) != None:
            return lodFlags[flag]
    return -1



def get_3dPosList_from_str(text):
    coor_list = [float(x) for x in text.split()]
    coor_list = [list(x) for x in zip(coor_list[0::3], coor_list[1::3], coor_list[2::3])]  # creating 2d coordinate array from 1d array
    return coor_list