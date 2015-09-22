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
from sklearn.datasets.base import Bunch
from sklearn.datasets import fetch_species_distributions
from sklearn.datasets.species_distributions import construct_grids
from sklearn import svm, metrics

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

class SupportVectorMachine(GeoAlgorithm):
    
    SPECIES = 'SPECIES'
    ENV = 'ENV'

    OUT_PRED = 'OUT_PRED'
    OUT_PRED_RES = 'OUT_PRED_RES'
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"svm.png")
    
    def defineCharacteristics(self):
        self.name = 'Support Vector Machine (SVM)'
        self.cmdName = 'svm'
        self.group = 'Species Distribution Modelling'

        self.addParameter(ParameterVector(self.SPECIES, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT,ParameterTable],False)) # Allow point
        self.addParameter(ParameterMultipleInput(self.ENV,'Environmental layers',ParameterMultipleInput.TYPE_RASTER, False))

        self.addOutput(OutputRaster(self.OUT_PRED,'Output Prediction'))
        self.addOutput(OutputTable(self.OUT_PRED_RES,'Stats'))

    def processAlgorithm(self, progress):
        # Set up the data as sklearn bunch (basically just a dictionary with specific attributes)
        data = Bunch()

        # Vector layer
        vector = self.getParameterValue(self.SPECIES)
        v = Processing.getObject(vector)
        v_crs = v.crs()
        
        # Environmental layers
        envlayers = self.getParameterValue(self.ENV)        
        if func.unificationNecessary(envlayers.split(";")):
            raise GeoAlgorithmExecutionException("All input environmental layers need to have the same resolution and extent. Use the Unify tool beforehand")
            #TODO: Enable option to do this automatically

        progress.setConsoleInfo("Loading Coverage Data")                

        # Check Projection and Cellsize
        for lay in envlayers.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            if r.crs() != v_crs:
                raise GeoAlgorithmExecutionException("All input layers need to have the same projection")
            if round(r.rasterUnitsPerPixelX()) != round(r.rasterUnitsPerPixelY()):
                raise GeoAlgorithmExecutionException("Grid Cell size values are not equal. Please be sure that grid cells are squares.")            

        # Set coverage parameters
        r = Processing.getObject(envlayers.split(";")[0]) # QgsRasterLayer object
        ex = r.extent()
        data["grid_size"] = r.rasterUnitsPerPixelX()        
        data["Nx"] = r.width()
        data["Ny"] = r.height()        
        data["x_left_lower_corner"] = ex.xMinimum()
        data["y_left_lower_corner"] = ex.yMinimum()

        # Load in Coverage values
        coverage = []
        for lay in envlayers.split(";"):
            raster = gdal.Open(str(lay))
            if raster.RasterCount > 1:
                progress.setConsoleInfo("Warning: Multiple bands for layer detected. Using only first band.")                
            array = raster.GetRasterBand(1).ReadAsArray()
            NA = raster.GetRasterBand(1).GetNoDataValue()
            if NA == None:
                raise GeoAlgorithmExecutionException("Warning: Raster layer has no no-data value. Please specify a no-data value for this dataset.")                
            else:
                array[array==NA] = -9999 # Replace nodata-values of array with -9999            
            coverage.append(array)    
        data["coverages"] = numpy.array( coverage ) # Load all the coverage values into the bunch

        # Setup parameters for output prediction
        a = gdal.Open(envlayers.split(";")[0])
        columns = a.RasterXSize
        rows = a.RasterYSize
        driver = a.GetDriver()            
        NA = -9999
        gt = a.GetGeoTransform()
        proj = a.GetProjection()
        output = self.getOutputValue(self.OUT_PRED)        


        # Set up the data grid
        xgrid, ygrid = construct_grids(data)
        
        # The grid in x,y coordinates           
        X, Y = numpy.meshgrid(xgrid, ygrid[::-1])

        # background points (grid coordinates) for evaluation        
        numpy.random.seed(100)
        background_points = numpy.c_[numpy.random.randint(low=0, high=data.Ny,
                                                    size=10000),
                                numpy.random.randint(low=0, high=data.Nx,
                                                    size=10000)].T

        # We'll make use of the fact that coverages[6] has measurements at all
        # land points.  This will help us decide between land and water.
        # FIXME: Assuming that all predictors have a similar distribution. Might be violated
        land_reference = data.coverages[0]

        progress.setConsoleInfo("Loading Occurence Data and coverage")                
        # Creating response
        train = []
        for feature in v.getFeatures():
            geom = feature.geometry().asPoint()
            mx = geom.x()
            my = geom.y()
            train.append((mx,my)) 
        data["train"] = numpy.array(train) # Add to bunch as training dataset

        # create species bunch      
        sp_Bunch = Bunch(name="Species")
        points = dict(train=data.train)
        for label, pts in points.iteritems():
            #determine coverage values for each of the training & testing points
            ix = numpy.searchsorted(xgrid, pts[0])
            iy = numpy.searchsorted(ygrid, pts[1])  
            bunch['cov_%s' % label] = data.coverages[:, -iy, ix].T

        progress.setConsoleInfo("Finished loading coverage data of environmental layers")                 
                                
        # Starting modelling
        progress.setConsoleInfo("Finished preparing the data for the analysis")                
        progress.setConsoleInfo("----")  
        progress.setConsoleInfo("Starting Modelling with support of sklearn")                
                      
        # Standardize features
        #TODO: Enable different or no Standardization methods
        mean = sp_Bunch.cov.mean(axis=0)
        std = sp_Bunch.cov.std(axis=0)
        train_cover_std = (sp_Bunch.cov - mean) / std

        # Fit OneClassSVM
        progress.setConsoleInfo("Fitting Support Vector Machine") 
        # TODO: Allow the user to vary the input                
        clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.5)
        clf.fit(train_cover_std)
        progress.setConsoleInfo("Fitting done") 
                
        # Predict species distribution using the training data
        Z = numpy.ones((data.Ny, data.Nx), dtype=numpy.float64)

        # We'll predict only for the land points.
        idx = numpy.where(land_reference > -9999)
        coverages_land = data.coverages[:, idx[0], idx[1]].T

        pred = clf.decision_function((coverages_land - mean) / std)[:, 0]
        Z *= pred.min()
        Z[idx[0], idx[1]] = pred

        levels = numpy.linspace(Z.min(), Z.max(), 25)
        Z[land_reference == -9999] = -9999

        result = Z # save the final results scores 
        
        # Compute AUC w.r.t. background points
        pred_background = Z[background_points[0], background_points[1]]
        pred_test = clf.decision_function((species.cov_test - mean)
                                          / std)[:, 0]
        scores = numpy.r_[pred_test, pred_background]
        y = numpy.r_[numpy.ones(pred_test.shape), numpy.zeros(pred_background.shape)]
        fpr, tpr, thresholds = metrics.roc_curve(y, scores)
        roc_auc = metrics.auc(fpr, tpr) #  Area under the ROC curve
        # TODO: Evaluate the availability of other metrics to compute on (average mean error, etc.. )
        
        # Create Output Prediction File
        output = self.getOutputValue(self.OUT_PRED_RES)
        titles =  ['AUC']
        res_pred = [roc_auc]
        # Save Output
        func.saveToCSV(res_pred, titles, output)

        # Create Output for resulting prediction        
        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Output creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            

        data_type = result.dtype        
        try:
            outData = driver.Create(output, columns, rows, 1, data_type)
        except Exception, e:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Output file could not be created!")
        
        band = outData.GetRasterBand(1)
        band.WriteArray(result)
        band.FlushCache()
        if nodata is not None:
            band.SetNoDataValue( nodata )
        band = None
    
        # Copy original geotransform and projection     
        outData.SetGeoTransform(gt) 
        outData.SetProjection(proj)

        outData = None # Close writing

    def helpFile(self):
        return os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
