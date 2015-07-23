# -*- coding: utf-8 -*-
"""
/***************************************************************************
 QSDM
        Species distribution modelling support for the QGIS Processing toolbox
                        -------------------
        begin                : 2014-03-31
        copyright            : (C) 2014 by Martin Jung
        email                : martinjung-at-zoho.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
__author__ = 'Martin Jung'
__date__ = 'April 2014'
__copyright__ = '(C) 2014, Martin Jung'
__revision__ = '$Format:%H$' # This will get replaced with a git SHA1 when you do a git archive

# Import PyQT bindings
from PyQt4.QtCore import *
from PyQt4.QtGui import *

# Import QGIS analysis tools
from qgis.core import *
from qgis.gui import *
from qgis.analysis import *
from qgis.utils import *

import os,sys,csv,string,math,operator,subprocess,tempfile,inspect

# Try to import functions from osgeo
try:
    from osgeo import gdal
except ImportError:
    import gdal
try:
    from osgeo import ogr, osr
except ImportError:
    import ogr,osr

import numpy

# Importing Sklean tools
import sklearn
from sklearn.linear_model import LogisticRegression, LinearRegression


from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.ProcessingLog import ProcessingLog
from processing.core.Processing import Processing
from processing.core.ProcessingConfig import ProcessingConfig
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
from processing.tools.system import *
# Also import settings
from qsdm_settings import qsdm_settings
# Import helperfunctions
import helperfunctions as func

# Import Processing algorithm parameters
try:
    from processing.parameters.ParameterTable import ParameterTable
    from processing.parameters.ParameterMultipleInput import ParameterMultipleInput
    from processing.parameters.ParameterRaster import ParameterRaster
    from processing.parameters.ParameterNumber import ParameterNumber
    from processing.parameters.ParameterSelection import ParameterSelection
    from processing.parameters.ParameterTableField import ParameterTableField
    from processing.parameters.ParameterExtent import ParameterExtent
    from processing.parameters.ParameterFixedTable import ParameterFixedTable
    from processing.parameters.ParameterVector import ParameterVector
    from processing.parameters.ParameterBoolean import ParameterBoolean
    from processing.parameters.ParameterFactory import ParameterFactory
    from processing.parameters.ParameterString import ParameterString
    
    from processing.outputs.OutputFactory import OutputFactory
    from processing.outputs.OutputTable import OutputTable
    from processing.outputs.OutputVector import OutputVector
    from processing.outputs.OutputRaster import OutputRaster
    from processing.outputs.Output import Output
except ImportError:
    from processing.core.parameters import ParameterTable
    from processing.core.parameters import ParameterMultipleInput
    from processing.core.parameters import ParameterRaster
    from processing.core.parameters import ParameterNumber
    from processing.core.parameters import ParameterSelection
    from processing.core.parameters import ParameterTableField
    from processing.core.parameters import ParameterExtent
    from processing.core.parameters import ParameterFixedTable
    from processing.core.parameters import ParameterVector
    from processing.core.parameters import ParameterBoolean
    #from processing.core.parameters import ParameterFactory
    from processing.core.parameters import ParameterString
    
    #from processing.core.outputs import OutputFactory
    from processing.core.outputs import OutputTable
    from processing.core.outputs import OutputVector
    from processing.core.outputs import OutputRaster
    from processing.core.outputs import Output  

class LogisticRegression(GeoAlgorithm):
    
    SPECIES = 'SPECIES'
    ENV = 'ENV'

    OUT_PRED = 'OUT_PRED'
    OUT_PRED_RES = 'OUT_PRED_RES'
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"logisticreg.png")
    
    def defineCharacteristics(self):
        self.name = 'Logistic Regression'
        self.cmdName = 'logisticreg'
        self.group = 'Species Distribution Modelling'

        self.addParameter(ParameterVector(self.SPECIES, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT,ParameterTable],False)) # Allow point
        self.addParameter(ParameterMultipleInput(self.ENV,'Environmental layers',ParameterMultipleInput.TYPE_RASTER, False))

        self.addOutput(OutputRaster(self.OUT_PRED,'Output Prediction'   ))
        self.addOutput(OutputTable(self.OUT_PRED_RES,'Stats'))

    def processAlgorithm(self, progress):
        # Vector layer
        vector = self.getParameterValue(self.SPECIES)
        v = Processing.getObject(vector)
        v_crs = v.crs()
        
        # Environmental layers
        envlayers = self.getParameterValue(self.ENV)        
        pred = []
        shape = None
        for lay in envlayers.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            if r.crs() != v_crs:
                raise GeoAlgorithmExecutionException("All input layers need to have the same projection")
            raster = gdal.Open(str(lay))
            gt = raster.GetGeoTransform()
            array = raster.GetRasterBand(1).ReadAsArray()
            if shape == None:
                shape = array.shape
            elif array.shape != shape:
                raise GeoAlgorithmExecutionException("All input environmental layers need to have the same resolution")      
            layers.append(array)

        # Creating response
        work_array = numpy.zeros_like(layers[0]) # Make a Presence/Absence of input array
        for feature in v.getFeatures():
            geom = feature.geometry().asPoint()
            mx = geom.x()
            my = geom.y()
            pp = func.world2Pixel(gt, mx,my)
            x = round(pp[0])
            y = round(pp[1])
            work_array[y,x] = 1
            
        response = work_array

        # The logistic regression
        logit = LogisticRegression(C=1.0,class_weight = 'auto',fit_intercept=True)
        logit.fit(train,response)


    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None
