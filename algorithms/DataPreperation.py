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
from processing.outputs.OutputFactory import OutputFactory
from processing.outputs.OutputTable import OutputTable
from processing.outputs.OutputVector import OutputVector
from processing.outputs.OutputRaster import OutputRaster
    

class CreateRichnessGrid(GeoAlgorithm):
    
    VECTOR = "VECTOR"
    SPEC_COL = "SPEC_COL"
    EXTENT = "EXTENT"
    GRAIN_SIZE = "GRAIN_SIZE"
    
    GRID = "GRID"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"default.png")
    
    def defineCharacteristics(self):
        self.name = 'Create Species Richness grid'
        self.cmdName = 'richnessgrid'
        self.group = 'Data Manipulation'
        
        self.addParameter(ParameterVector(self.VECTOR, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT],False)) # Allow point and later polygons ,ParameterVector.VECTOR_TYPE_POLYGON
        self.addParameter(ParameterTableField(CreateRichnessGrid.SPEC_COL, "Species Name", CreateRichnessGrid.VECTOR))
        self.addParameter(ParameterNumber(self.GRAIN_SIZE, 'Grain size',1,None,100))
        self.addParameter(ParameterExtent(self.EXTENT, 'Extent of the new Grid'))
        
        self.addOutput(OutputRaster(self.GRID,'Richness Grid'))
    
    def processAlgorithm(self, progress):
        # Do the stuff
        vector = self.getParameterValue(self.VECTOR)
        v = Processing.getObject(vector)
        scl = self.getParameterValue(self.SPEC_COL)
        ext = self.getParameterValue(self.EXTENT)
        try:
            ext = string.split(ext,",") # split 
        except AttributeError: # Extent was empty, raise error
            raise GeoAlgorithmExecutionException("Please set an extent for the generated raster")        
        cs  =  self.getParameterValue(self.GRAIN_SIZE)
        output = self.getOutputValue(self.GRID)

        # Create output layer
        xmin = float(ext[0])
        xmax = float(ext[1])
        ymin = float(ext[2])
        ymax = float(ext[3])
        gt = (xmin,cs,0.0,ymax,0.0,-cs)
        nodata = -9999
        
        cols = int( abs( (xmax-xmin)/gt[1] ) )
        rows = int( abs( (ymax-ymin)/gt[5] ) )
        fin_array = numpy.zeros((rows,cols)) # Create empty grid

        #if vector is a point do the following, else calculate for overlapping range sizes
        if v.geometryType() == QGis.Point:
            progress.setConsoleInfo("Using the point layers to calculate Species richness for resulting grid.")
            progress.setConsoleInfo("---")
            # Get the number of species
            noSpecies = func.getUniqueAttributeList( v, scl)
            progress.setConsoleInfo("Processing %s number of different species" % (str(len(noSpecies))) )

            ds = ogr.Open(vector)
            name = ds.GetLayer().GetName()
            proj = ds.GetLayer().GetSpatialRef()
            n = ds.GetLayer().GetFeatureCount()

            k = 1
            for spec in noSpecies:
                # Make a copy of the final_array
                work_array = numpy.zeros_like(fin_array)
                # Vector layer subsetting to the specific species
                layers = ds.ExecuteSQL("SELECT * FROM %s WHERE %s = '%s'" % (name, scl, spec) )
                progress.setConsoleInfo("Gridding %s individual points of species %s " % (str(layers.GetFeatureCount()), spec ))                
                func.updateProcessing(progress,k,n )
                for i in range(0,layers.GetFeatureCount()):
                    f = layers.GetFeature(i)
                    geom = f.GetGeometryRef()
                    mx,my= geom.GetX(), geom.GetY()  #coord in map units                
                    pp = func.world2Pixel(gt, mx,my)                
                    x = round(pp[0])
                    y = round(pp[1])

                    if x < 0 or y < 0 or x >= work_array.shape[1] or y >= work_array.shape[0]:
                        progress.setConsoleInfo("Point %s outside given exent" % (str( f.GetFID() )) )
                    else:            
                        # Check if species was already added to grid cell
                        test = work_array[y,x]
                        if test != 1:
                            work_array[y,x] = 1
                k += 1
                # Add the working arrays values to 
                fin_array = numpy.add(work_array,fin_array)

        elif v.geometryType() == QGis.Polygon:
            progress.setConsoleInfo("Using the range size polygons to calculate Species richness for resulting grid.")

            # rasterization        
        if numpy.count_nonzero(fin_array) == 0:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"No values were rasterized. Check GeometryType and Vector Projection.")
            
            
        # Create output raster
        func.createRaster(output,cols,rows,fin_array,nodata,gt,proj,'GTiff')
        # And free up memory
        del(ds,layers)
        

class DataTransformationSimple(GeoAlgorithm):
    
    RASTER = "RASTER"  
    TRANSFORM = "TRANSFORM"
    tsel = ["Log10","Log-e","Log-e+1","Square-Root","* -1","Exp-e","ArcSin","Polynom n","Standardize: (X - Mean)/SD","Normalize: [0,1]"]
    POLYNOM = "POLYNOM"
    
    OUT = "OUT"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"transform.png")
    
    def defineCharacteristics(self):
        self.name = 'Data transformation (Simple)'
        self.cmdName = 'transformationSimple'
        self.group = 'Data Manipulation'
        
        self.addParameter(ParameterRaster(self.RASTER, "Input Grid", False))
        self.addParameter(ParameterSelection(self.TRANSFORM, "Transformation method", self.tsel, 0))
        self.addParameter(ParameterNumber(self.POLYNOM, "Polynominal Expansion on n=", 2, None, 1))
        
        self.addOutput(OutputRaster(self.OUT,'Transformed Raster'))
    
    def processAlgorithm(self, progress):
        # Do the thing, Julie
        inputFilename = self.getParameterValue(self.RASTER)
        what = self.tsel[self.getParameterValue(self.TRANSFORM)]
        polyn = self.getParameterValue(self.POLYNOM)
        
        output = self.getOutputValue(self.OUT)
        # Starting transformation
        progress.setConsoleInfo("Transforming input raster layer with %s " % (str( what )))
        func.updateProcessing(progress,1,3)
 
        raster = gdal.Open(str(inputFilename))
        columns =  raster.RasterXSize
        rows = raster.RasterYSize
        driver = raster.GetDriver()
        if(raster.RasterCount==1):
            band = raster.GetRasterBand(1)
            data_type = band.DataType
            nodata = band.GetNoDataValue()          
            # Raise Notice that no no-data value has been selected
            progress.setConsoleInfo("!!! - Found no nodata-value for band %s . Transforming all values - !!!" % (str(1)) )
            try:
                array = band.ReadAsArray() 
            except Exception, e:
                ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Could not transform the raster to an array.")
            work = array[array!=nodata]
            res = self.transformArray(work,what,polyn) # Transformation
            result = numpy.copy(array)
            result[result!=nodata] = res  # Replace new values with the one sampled
            func.updateProcessing(progress,2,3)
        else:
            result = []
            for i in range(0,raster.RasterCount):
                band = raster.GetRasterBand(i)
                data_type = band.DataType
                nodata = band.GetNoDataValue()          
                # Raise Notice that no no-data value has been selected
                progress.setConsoleInfo("!!! - Found no nodata-value for band %s . Transforming all values - !!!" % (str(i+1)) )
                try:
                    array = band.ReadAsArray() 
                except Exception, e:
                    ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Could not transform the raster to an array.")
                work = array[array!=nodata]
                tran = self.transformArray(work,what,polyn)
                r = numpy.copy(array)
                r[r!=nodata] = tran  # Replace new values with the one sampled
                result.append(r)
            
        # Create Output
        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            
        
        try:
            outData = driver.Create(output, columns, rows, raster.RasterCount, data_type)
        except Exception, e:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Output file could not be created!")
        
        for i in range(0,raster.RasterCount):
            band = outData.GetRasterBand(i+1)
            if raster.RasterCount == 1:
                band.WriteArray(result)
            else:
                band.WriteArray(result[i]) 
            band.FlushCache()
            na = raster.GetRasterBand(i+1).GetNoDataValue()
            if na is not None:
                band.SetNoDataValue( nodata )
            band = None
    
        func.updateProcessing(progress,3,3)    
        # Copy original geotransform and projection     
        outData.SetGeoTransform(raster.GetGeoTransform()) 
        outData.SetProjection(raster.GetProjection())
        
        outData = None # Close writing

        # Last check
        if os.path.exists(output) == False:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Transformation could not be generated.")
            
    def transformArray(self, array, what,polyn=2):
        if what == "Log10":
            res = numpy.log10(array)
            return res
        elif what == "Log-e":
            res = numpy.log(array)
            return res
        elif what == "Log-e+1":
            res = numpy.log1p(array)
            return res
        elif what == "Square-Root":
            res = numpy.sqrt(array)
            return res
        elif what == "* -1":
            res = numpy.multiply(array,-1)
            return res
        elif what == "Exp-e":
            res = numpy.exp(array)
            return res
        elif what == "ArcSin":
            res = numpy.arcsin(array)
            return res
        elif what == "Polynom n":
            res = numpy.polynomial.polynomial.polymul(array,polyn)
            return res
        elif what == "Standardize: (X - Mean)/SD":
            res = ( array - numpy.nanmean(array,axis=0)) / numpy.nanstd(array,axis=0,ddof=1)
            return res
        elif what == "Normalize: [0,1]":
            res = (array - numpy.nanmin(array,axis=0)) / (numpy.nanmax(array,axis=0) - numpy.nanmin(array,axis=0))
            return res

# Mostly inspired by this brilliant posting 
# http://kldavenport.com/mahalanobis-distance-and-outliers/
class VectorOutlierSelection(GeoAlgorithm):
    """
    Detects Outliers in Point data based on their Mahalanobis distances
    """    
    VECTOR = "VECTOR"
#    MULTIPLIER = "MULTIPLIER"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"outlier.png")
    
    def defineCharacteristics(self):
        self.name = 'Select Outlier Points'
        self.cmdName = 'outlierselect'
        self.group = 'Data Manipulation'
        
        self.addParameter(ParameterVector(self.VECTOR, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT],False))
#        self.addParameter(ParameterNumber(self.MULTIPLIER, 'Multiplier',1,None,2.0))
    
    def processAlgorithm(self, progress):
        # Do the stuff
        self.progress = progress
        vector = self.getParameterValue(self.VECTOR)
        v = Processing.getObject(vector)
#        mp  =  self.getParameterValue(self.MULTIPLIER)

        # Get List of coordinates
        progress.setConsoleInfo("Get Input coordinates...")
        func.updateProcessing(progress,1,4)    

        x = y = id = []
        ds = ogr.Open(vector)
        lay = ds.GetLayer()
        for i in range(0,lay.GetFeatureCount()):
            f = lay.GetFeature(i)
            geom = f.GetGeometryRef()
            x.append(geom.GetX())
            y.append(geom.GetY())
            id.append(f.GetFID())

        if len(x) == 0 or len(y) == 0:
            raise GeoAlgorithmExecutionException("Coordinates of given point layer could not be extracted")      

        # Build Mahalanobis Distances
        progress.setConsoleInfo("Build Mahalanobis Distances...")
        func.updateProcessing(progress,2,4)    

        md = self.MahalanobisDist(x,y)
        
        # Identify outliers and build new values
        progress.setConsoleInfo("Identify outliers...")
        func.updateProcessing(progress,3,4)    

        # Get 3 greatest values indices
        outliers = (-numpy.array(md)).argsort()[:3] 

#         threshold = numpy.mean(md) * mp # adjust 1.5 accordingly 
#         outliers = []
#         for i in range(len(md)):
#             if md[i] >= threshold:
#                 outliers.append(i) # position of removed pair
        
        if len(outliers) >= (len(id) / 4):
            raise GeoAlgorithmExecutionException("Too many outliers. Try to increase the multiplier.")
        
        # Get ids of outliers and select those in the input vectorlayer
        progress.setConsoleInfo("Select %s outliers in the input vectorlayer" % (str(len(outliers))))
        func.updateProcessing(progress,4,4)
         
        out = []
        for o in outliers:
            i = int( id[o] )
            v.select(i)
        
    def MahalanobisDist(self,x, y):
        # Estimate a covariance matrix for (x,y)
        covariance_xy = numpy.cov(x,y, rowvar=0)
        try:
            inv_covariance_xy = numpy.linalg.inv(covariance_xy)
        except numpy.linalg.LinAlgError:
            # There is no linear inverse matrix for given points
            self.progress.setConsoleInfo("Singular matrix. Calculating (Moore-Penrose) pseudo-inverse matrix instead.")
            inv_covariance_xy = numpy.linalg.pinv(covariance_xy)
            #raise GeoAlgorithmExecutionException("Singular non-invertable covariance matrix. Looking for solutions for this.")      

        # Center each value by the mean
        xy_mean = numpy.mean(x),numpy.mean(y)
        x_diff = numpy.array([x_i - xy_mean[0] for x_i in x])
        y_diff = numpy.array([y_i - xy_mean[1] for y_i in y])
        diff_xy = numpy.transpose([x_diff, y_diff])
        
        # Formula for MahalanobisDist
        md = []
        for i in range(len(diff_xy)):
            md.append(numpy.sqrt(numpy.dot(numpy.dot(numpy.transpose(diff_xy[i]),inv_covariance_xy),diff_xy[i])))
        return md
    
