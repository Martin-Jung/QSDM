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
#import numpy.lib.recfunctions as rfn # numpy.recfunctions

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
    from processing.outputs.OutputFactory import OutputFactory
    from processing.outputs.OutputTable import OutputTable
    from processing.outputs.OutputVector import OutputVector
    from processing.outputs.OutputRaster import OutputRaster
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
    #from processing.core.outputs import OutputFactory
    from processing.core.outputs import OutputTable
    from processing.core.outputs import OutputVector
    from processing.core.outputs import OutputRaster


class NicheOverlapStatistics(GeoAlgorithm):
    
    ENV = 'ENV'
    METRIC = 'METRIC'
    m = ['Schoener\'s D','Warren\'s I']
    
    RESULTS = "RESULTS"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"nicheoverlap.png")
    
    def defineCharacteristics(self):
        self.name = 'Calculate Niche Overlap Statistics'
        self.cmdName = 'nicheoverlap'
        self.group = 'Data Analysis'
        
        self.addParameter(ParameterMultipleInput(self.ENV,'Prediction layers',ParameterMultipleInput.TYPE_RASTER, False))
        self.addParameter(ParameterSelection(self.METRIC, "Metric", self.m, 0))
        
        self.addOutput(OutputTable(self.RESULTS,'Result output'))
    
    def processAlgorithm(self, progress):
        # Do the stuff
        metric = self.m[self.getParameterValue(self.METRIC)]
        envlayers = self.getParameterValue(self.ENV)        
        names = []
        env = []  
        shape = None
        # Load all arrays into an dictonary
        for lay in envlayers.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            name = str( r.name() )
            names.append(name)
            a = gdal.Open( lay )
            array = a.GetRasterBand(1).ReadAsArray()
            # Prepare by removing all no-data values from array
            NA = a.GetRasterBand(1).GetNoDataValue()
            if NA != None:
                array[array==NA] = 0
            env.append(array)            
            if shape == None: # Fast check if array are unequal
                shape = array.shape
            else:
                if shape != array.shape:
                    raise GeoAlgorithmExecutionException("Input layers need to have the same extent and cellsize.") 
            a = r = None
        if len( env ) == 1:
            raise GeoAlgorithmExecutionException("You need at least two layers to calculate overlap statistics.")             
        
        progress.setConsoleInfo("Loaded %s arrays for calculation" % ( str( len( names) ) ) )
        func.updateProcessing(progress,1,3)
        
        results = []
        func.updateProcessing(progress,2,3)
        if len(env) > 2:
        # Iterative Calculation of the metrics
            for j in range(0,len(env)):
                for k in range(j+1,len(env)):
                    progress.setConsoleInfo("Calculating Overlap of layers %s with %s" % (names[j],names[k]) )
                    r = self.Overlap(metric,env[j],env[k])
                    res = (names[j],names[k],metric,r)
                    results.append(res)
        else:
            # Only two input layers
            r = self.Overlap(metric,env[0],env[1])
            res = (names[0],names[1],metric,r)
            results.append(res)
            
        progress.setConsoleInfo("Saving results")
        func.updateProcessing(progress,3,3)
        
        output = self.getOutputValue(self.RESULTS)
        titles =  ['Layer1','Layer2','Metric','Overlap']
        # Save Output
        func.saveToCSV(results, titles, output )
        
    def Overlap(self,m,array1,array2):
        a1_sum = numpy.sum(array1)
        a2_sum = numpy.sum(array2)
        if m == 'Schoener\'s D':
            # Following Warren 2008, after Schoener 1968
            D = numpy.sum( numpy.abs( numpy.double(array1)/a1_sum - numpy.double(array2)/a2_sum ) )
            res = 1 - 0.5 * D
            return res            
        elif m == 'Warren\'s I':
            # Following Warren 2008 + Erratum using Hellingers Distances
            # Calculate hellinger distance after van der Vaart 1998
            H = numpy.power(numpy.sqrt( numpy.double(array1)/a1_sum) - numpy.sqrt( numpy.double(array2)/a2_sum),2)
            # Following Warrens erratum -> 10.1111/j.1558-5646.2010.01204.x
            res = 1 - ( numpy.sum(H)  * 0.5 )
            return res
        else:
            return None

    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None

# MESS calculation after Elith 2010
# Code adapted from Jean-Pierre Rossi, Robert Hijmans, Paulo van Breugel 
class MESS(GeoAlgorithm):

    SPECIES = 'SPECIES'
    ENV = 'ENV'    
    RESULTS = "RESULT"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"default.png")
    
    def defineCharacteristics(self):
        self.name = 'Multivariate Environmental Similarity Surfaces'
        self.cmdName = 'mess'
        self.group = 'Data Analysis'

        self.addParameter(ParameterVector(self.SPECIES, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT,ParameterTable],False))        
        self.addParameter(ParameterMultipleInput(self.ENV,'Environmental layers',ParameterMultipleInput.TYPE_RASTER, False))
        
        self.addOutput(OutputRaster(self.RESULTS,'MESS'))
    
    def processAlgorithm(self, progress):
        # Do the stuff
        point = self.getParameterValue(self.SPECIES)
        v = Processing.getObject(point) 
        ref_crs = v.crs() # For Comparison with the raster layers
        x = y = v_id = []
        # Iter through and get coordinates and id
        iter = v.getFeatures()
        for feature in iter:
            v_id.append( feature.id() )
            geom = feature.geometry().asPoint()
            x.append(geom.x())
            y.append(geom.y())
        
        envlayers = self.getParameterValue(self.ENV)        
        names = []
        env = []  
        shape = None
        # Load all arrays into an dictonary
        for lay in envlayers.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            name = str( r.name() )
            if r.crs() != ref_crs:
                raise GeoAlgorithmExecutionException("Input point layer and all environmental layers need to have the same projection.")                 
            names.append(name)
            a = gdal.Open( lay )
            array = a.GetRasterBand(1).ReadAsArray()
            env.append(array)            
            if shape == None: # Fast check if array are unequal
                shape = array.shape
            else:
                if shape != array.shape:
                    raise GeoAlgorithmExecutionException("Input layers need to have the same extent and cellsize.") 
            a = r = None
        if len( env ) == 1:
            raise GeoAlgorithmExecutionException("You need more than two layers to calculate the MESS index.")             
        
        progress.setConsoleInfo("Loaded %s arrays for calculation" % ( str( len( names) ) ) )
        func.updateProcessing(progress,1,5)
        
        ## Start Calculating MESS
        # Get Point Values
        values = self.extractPointValues(envlayers,names,x,y)
        progress.setConsoleInfo("Extracted Point data values")
        func.updateProcessing(progress,2,5)

        # Flattening predictor values
        r = numpy.array(env)
        ref = r.reshape((r.shape[0], -1))         
        
        # Flattening extracted point values
        v = numpy.array(values)
        val = v.reshape((v.shape[0], -1))  

        func.updateProcessing(progress,2,4)
        # Iterative Calculation of the metric
        out = []
#       for i in range(0,len(env)):
#            numpy.apply_along_axis(self.getMESS,1,ref,ref[,i])
        result = None

        #         mess<-function(X,V,full=TRUE)

        # for (i in 1:(dim(E)[2])) {
        #     e<-data.frame(E[,i]) ; v<-V[,i]
        #     r_mess[[i]][]<-apply(X=e, MARGIN=1, FUN=messi, v=v)
        #     }
        # rmess<-r_mess[[1]]
        # E<-extract(x=r_mess,y=1:ncell(r_mess[[1]]))
        # rmess[]<-apply(X=E, MARGIN=1, FUN=min)
        # if(full==TRUE) {
        #     out     layerNames(out)<-c(layerNames(X),"mess")
        # }
        # if(full==FALSE) out return(out)
        # }
        
        
        progress.setConsoleInfo("Saving results")
        func.updateProcessing(progress,5,5)
        
        # Create Output
        a = gdal.Open(envlayers.split(";")[0])
        columns = a.RasterXSize
        rows = a.RasterYSize
        driver = a.GetDriver()            
        NA = a.GetRasterBand(1).GetNoDataValue()
        data_type = a.GetRasterBand(1).DataType
        gt = a.GetGeoTransform()
        proj = a.GetProjection()
        output = self.getOutputValue(self.RESULTS)        
        
        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            
        
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
    
    def extractPointValues(self,envlayers,x,y):
        ''' 
        Get the point values from all input rasters
        Returns an array with all values for the rasters
        '''
        res = []
        c = 0
        for lay in envlayers.split(";"):
            raster = gdal.Open(lay)
            gt = raster.GetGeoTransform() # Get geotransform
            array = raster.GetRasterBand(1).ReadAsArray()
            raster = None # close gdal
            vals = []
            for i in range(0, len(x)):
                pp = func.world2Pixel(gt, x[i],y[i])
                mx = int(pp[0])
                my = int(pp[1])
                vals.append( array[y, x] )
            res.append(vals)
            c += 1
        return res

    #  Plain conversion to Python from Jean-Pierre Rossi R code
    #  Could probably be optimized
    def getMESS(self,p,v):
        niinf = len( numpy.where(v<=p) )
        f =  100 * niinf / len(v)
        if (f == 0):
            simi = 100 * ( p - min( v ) ) / ( max( v ) - min( v ))
        if (0 < f and f <= 50):
            simi = 2 * f
        if (50 <= f and f < 100):
            simi = 2 * (100 - f)
        if (f == 100):
            simi = 100 * ( max( v ) - p) / ( max( v ) - min( v ))
        return simi

        # .messi2 <- function(p,v){
        # 	v <- na.omit(v)
        # 	f <- 100*findInterval(p, sort(v)) / length(v)
        # 	minv <- min(v)
        # 	maxv <- max(v)
        # 	ifelse(f == 0, 100*(p-minv)/(maxv-minv), 
        # 		ifelse(f <= 50, 2*f, 
        # 		ifelse(f < 100, 2*(100-f),
        # 			100*(maxv-p)/(maxv-minv)
        # 	)))
        # }


        
    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None

# Following the quantitative approach from Mesgaran et al. 2014
class NovelConditions(GeoAlgorithm):
    
    REFERENCE = 'REFERENCE'
    PROJECTION = 'PROJECTION'
    RESULT = "RESULT"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"novelcondition.png")
    
    def defineCharacteristics(self):
        self.name = 'Evaluate Novel Conditions'
        self.cmdName = 'novelcondition'
        self.group = 'Data Analysis'

        self.addParameter(ParameterMultipleInput(self.REFERENCE,'Reference layers',ParameterMultipleInput.TYPE_RASTER, False))
        self.addParameter(ParameterMultipleInput(self.PROJECTION,'Projection layers',ParameterMultipleInput.TYPE_RASTER, False))
        self.addOutput(OutputRaster(self.RESULT,'Novel Conditions Map'))
    
    def processAlgorithm(self, progress):
        # Do the stuff, julie
        reference = self.getParameterValue(self.REFERENCE)        
        projection = self.getParameterValue(self.PROJECTION)        
        ref_ar = []
        proj_ar = []  
        crs = None
        shape = None
        # Load all reference arrays into an dictonary
        for lay in reference.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            if crs == None:
                crs = r.crs()
            else:
                if r.crs() != crs:
                    raise GeoAlgorithmExecutionException("All input layers need to have the same projection.")                 
            a = gdal.Open( lay )
            array = a.GetRasterBand(1).ReadAsArray()
            ref_ar.append(array)
            if shape == None: # Fast check if array are unequal
                shape = array.shape
            else:
                if shape != array.shape:
                    raise GeoAlgorithmExecutionException("Input layers need to have the same extent and cellsize.") 
            a = r = None
        # Load all projection arrays into an dictonary
        for lay in projection.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            if crs == None:
                crs = r.crs()
            else:
                if r.crs() != crs:
                    raise GeoAlgorithmExecutionException("All input layers need to have the same projection.")                 
            a = gdal.Open( lay )
            array = a.GetRasterBand(1).ReadAsArray()
            proj_ar.append(array)
            if shape == None: # Fast check if array are unequal
                shape = array.shape
            else:
                if shape != array.shape:
                    raise GeoAlgorithmExecutionException("Input layers need to have the same extent and cellsize and nodata value.") 
            a = r = None
            
        
        if len( ref_ar ) != len( proj_ar ):
            raise GeoAlgorithmExecutionException("Need the same number of reference layers as projection layers.")            

        progress.setConsoleInfo("Successfully loaded layers for calculation" )
        func.updateProcessing(progress,1,5)

        # Take Values for output from the first reference layer
        a = gdal.Open( reference.split(";")[0] )
        columns = a.RasterXSize
        rows = a.RasterYSize
        driver = a.GetDriver()            
        data_type = a.GetRasterBand(1).DataType
        gt = a.GetGeoTransform()
        proj = a.GetProjection()
        nodata = a.GetRasterBand(1).GetNoDataValue()
        output = self.getOutputValue(self.RESULT)        
        a = None
                

        # test
        a = "/home/martin/.qgis2/python/plugins/QSDM/sampledata/environment/cur/t_mean.asc"
        a = "/home/martin/.qgis2/python/plugins/QSDM/sampledata/environment/cur/elev_mean.asc"
        r = gdal.Open(a)
        ref_ar = r.GetRasterBand(1).ReadAsArray()
        proj_ar = r.GetRasterBand(1).ReadAsArray()
        nodata = -9999
        
        ref = ref_ar.ravel()
        proj = proj_ar.ravel()
        # Reshaping to a single shape
        progress.setConsoleInfo("Flattening Input layers to a single shape" )
        func.updateProcessing(progress,2,5)        
        r = numpy.array(ref_ar)# Reference
        ref = r.reshape((r.shape[0], -1))         
        r = numpy.array(proj_ar)# Projection
        proj = r.reshape((r.shape[0], -1))         
        
        # Calculating Average and Covariance
        progress.setConsoleInfo("Calculating Average and Covariance" )        
        func.updateProcessing(progress,3,5)
        ref = ref.astype(float)        
        ref[ref == nodata] = numpy.nan # Replace values of nodata with NaN
        
        ref_avg = numpy.apply_along_axis(numpy.nanmean,1,ref)
            
        m = numpy.ma.make_mask((numpy.isnan(ref)))
        mask = numpy.ma.MaskedArray(ref,m) # use a mask to kickout the nan values
        ref_cov = numpy.ma.cov(ref,mask,rowvar=0,allow_masked=True)

        numpy.ma.reshape(ref_cov, (len(ref_avg),-1))
        
        try:    
            inv_covariance_xy = numpy.linalg.inv(ref_cov.data)       
        except numpy.linalg.LinAlgError:
            # There is no linear inverse matrix for given points
            self.progress.setConsoleInfo("Singular matrix. Calculating (Moore-Penrose) pseudo-inverse matrix instead.")
            #inv_covariance_xy = numpy.linalg.pinv(ref_cov)
            raise GeoAlgorithmExecutionException("Singular non-invertable covariance matrix. Looking for solutions for this.")      


        # Calculate Mahalanobis distance ratios
        progress.setConsoleInfo("Calculating Mahalanobis ratios between raster cells")
        func.updateProcessing(progress,4,5)

        try:
            from scipy.spatial.distance import mahalanobis
            
            # Calculate for ref arrays
            mah.ref = mahalanobis(ref,ref_avg,inv_covariance_xy)
            mah.ref = numpy.exp2(mah.ref) # Calculate D^2 - squared
            
            # Calculate for proj arrays
            mah.proj = mahalanobis(proj,ref_avg,inv_covariance_xy)
            mah.proj = numpy.exp2(mah.proj) # Calculate D^2 - squared
            

            # Ratios
            mah.max = numpy.max( mah.ref[numpy.isfinite( mah.ref )] )
            result <- numpy.divide(mah.proj,mah.max )            
            
        except ImportError:
            self.progress.setConsoleInfo("Scipy not found. Calculating mahalanobis manually.")
            
            # Center each value by the mean
            x = ref_ar[0]
            xy_mean = ref_avg
            x_diff = numpy.array([x_i - xy_mean[0] for x_i in x])
            y_diff = numpy.array([y_i - xy_mean[1] for y_i in y])
            diff_xy = numpy.transpose([x_diff, y_diff])       
            # Formula for MahalanobisDist
            md = []
            for i in range(len(diff_xy)):
                md.append( numpy.sqrt(numpy.dot(numpy.dot(numpy.transpose(diff_xy[i]),inv_covariance_xy),diff_xy[i])) )
            return md
    
        progress.setConsoleInfo("Saving results")
        func.updateProcessing(progress,5,5)

        # Create Output
        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            
        
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
        
        # Load in QGIS
        func.rasterInQgis(output)

    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None

class RangeShifts(GeoAlgorithm):
    
    BASELINE = 'BASELINE'
    PREDICTION = 'PREDICTION'
    QUANT = 'QUANT'
    RESULTS = "RESULT"
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"rangeshift.png")
    
    def defineCharacteristics(self):
        self.name = 'Range Shifts'
        self.cmdName = 'rangeshift'
        self.group = 'Data Analysis'

        self.addParameter(ParameterRaster(self.BASELINE, "Baseline Grid", False))
        self.addParameter(ParameterRaster(self.PREDICTION, "Prediction Grid", False))        
        self.addParameter(ParameterNumber(self.QUANT, "Quantile distance from Median", 10, None, 1))
        
        self.addOutput(OutputRaster(self.RESULTS,'Range Shift Map'))
    
    def processAlgorithm(self, progress):
        # Do the stuff
        b = self.getParameterValue(self.BASELINE)        
        p = self.getParameterValue(self.PREDICTION)        
        d = self.getParameterValue(self.QUANT)
        
        # baseline
        r = Processing.getObject(b)
        baseline_name = str( r.name() )
        a = gdal.Open( b )
        baseline = a.GetRasterBand(1).ReadAsArray()
        # Prepare by removing all no-data values from array
        nodata = a.GetRasterBand(1).GetNoDataValue()
        if nodata == None:
            raise GeoAlgorithmExecutionException("Please classify both layers with a valid nodata-value.") 

        # Take Values for output from baseline raster
        columns = a.RasterXSize
        rows = a.RasterYSize
        driver = a.GetDriver()            
        data_type = a.GetRasterBand(1).DataType
        gt = a.GetGeoTransform()
        proj = a.GetProjection()
        output = self.getOutputValue(self.RESULTS)        
        a = None
        
        # Prediction
        r = Processing.getObject(p)
        prediction_name = str( r.name() )
        a = gdal.Open( b )
        prediction = a.GetRasterBand(1).ReadAsArray()
        # Prepare by removing all no-data values from array
        NA = a.GetRasterBand(1).GetNoDataValue()
        if NA == None:
            raise GeoAlgorithmExecutionException("Please classify both layers with a valid nodata-value.") 
        a = None
        
        # Compare shapes
        if baseline.shape != prediction.shape:
            raise GeoAlgorithmExecutionException("Input layers need to have the same extent and cellsize.") 
        
        progress.setConsoleInfo("Loaded layers for calculation" )
        func.updateProcessing(progress,1,4)
        
        progress.setConsoleInfo("Starting processing" )
        func.updateProcessing(progress,2,4)

        # Index values
        # 1  = Range contraction
        # 2  = No Change
        # 3  = Range Expansion
        med = numpy.median(baseline[baseline!=nodata])
        low = numpy.percentile(baseline[baseline!=nodata],50-d)
        high = numpy.percentile(baseline[baseline!=nodata],50+d)
        
        res = numpy.zeros_like(baseline)
        res[numpy.where(numpy.logical_and(prediction >= low,prediction <= high))] = 2
        res[prediction < low] = 1
        res[prediction > high] = 3
        res[baseline==nodata] = nodata # fill the rest with nodata
        
        result = numpy.copy(res)
        progress.setConsoleInfo("Saving results")
        func.updateProcessing(progress,3,4)

        # Create Output
        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            
        
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
        
        # Load in QGIS
        func.rasterInQgis(output)
        
        progress.setConsoleInfo("Styling")
        func.updateProcessing(progress,4,4)
        if os.path.exists(output) == False:
            raise GeoAlgorithmExecutionException("Change Map could not be generated.") 

        ## Styling 
        canvas = QgsMapCanvas()
        canvas.freeze(True)           
        lyr = func.getLayerByName("RESULT")
        lyr.setDrawingStyle("SingleBandPseudoColor")
        # The band of classLayer
        classLyrBnd = 1
        # Color list for ramp
        clrLst = [ QgsColorRampShader.ColorRampItem(1, QColor(255,51,51),"Niche contraction"),  # Contraction red
                    QgsColorRampShader.ColorRampItem(2,QColor(255,255,192),"No change"), # no change, yellowish
                    QgsColorRampShader.ColorRampItem(3, QColor(0,0,255),"Niche expansion") ] # expansion, blue
        #Create the shader
        lyrShdr = QgsRasterShader()
        #Create the color ramp function
        clrFnctn = QgsColorRampShader()
        clrFnctn.setColorRampType(QgsColorRampShader.EXACT)
        clrFnctn.setColorRampItemList(clrLst)
        #Set the raster shader function
        lyrShdr.setRasterShaderFunction(clrFnctn)
        #Create the renderer
        lyrRndr = QgsSingleBandPseudoColorRenderer(lyr.dataProvider(), classLyrBnd, lyrShdr)
        #Apply the renderer to classLayer
        lyr.setRenderer(lyrRndr)
        #refresh legend
        if hasattr(lyr, "setCacheImage"):
            lyr.setCacheImage(None)
        lyr.triggerRepaint()
        iface.legendInterface().refreshLayerSymbology(lyr)
        canvas.freeze(False)
        canvas.refresh()

    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None
            