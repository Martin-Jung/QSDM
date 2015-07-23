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

import os, subprocess, ntpath, platform
import numpy
try:
    from osgeo import gdal
except ImportError:
    import gdal

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


class Maxent(GeoAlgorithm):
    
    SPECIES = 'SPECIES'
    SPEC_COL = 'SPEC_COL'
    ENV = 'ENV'
    #ENV_PRED = 'ENV_PRED'
    PARAM = 'PARAM'

#    OUTPUT_DIR = 'OUTPUT_DIR'
#    OUT_PRED = 'OUT_PRED'
#    OUT_PRED_RES = 'OUT_PRED_RES'
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"maxent.png")
    
    def defineCharacteristics(self):
        self.name = 'Maximum Entropy Modelling'
        self.cmdName = 'maxent'
        self.group = 'Species Distribution Modelling'

        self.addParameter(ParameterVector(self.SPECIES, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT,ParameterTable],False)) # Allow point
        self.addParameter(ParameterTableField(Maxent.SPEC_COL, "Species Name", Maxent.SPECIES))
        self.addParameter(ParameterMultipleInput(self.ENV,'Environmental layers',ParameterMultipleInput.TYPE_RASTER, False))
        self.addParameter(ParameterTable(self.PARAM, '(Optional) MAXENT parameter file'))
        #self.addParameter(ParameterMultipleInput(self.ENV_PRED,'Projection of Environmental layers (Optional)',ParameterMultipleInput.TYPE_RASTER, True))

#        self.addOutput(Output(self.OUTPUT_DIR,'MAXENT output folder'   ))
#        self.addOutput(OutputTable(self.OUT_PRED_RES,'Prediction Accuracy Results'))
    
    def processAlgorithm(self, progress):
        ## Parameter preperation ##
        # Get the location to the maxent jar file
        maxent = qsdm_settings.maxent()
        if os.path.basename(maxent) == 'maxent.jar':
            pass
        else:
            maxent = os.path.join(qsdm_settings.maxent(), 'maxent.jar') # If the directory and not the file was chosen
        # Get location of java, available memory and output path
        if sys.platform == "win32" or "win64":
            ex = "java.exe"
        else:
            ex = "java"
        java = os.path.join(qsdm_settings.javaPath(),ex) # the path to java if basic execution fails
        mem = str( qsdm_settings.getMEM() )     # available memory for MAXENT
        work = qsdm_settings.workPath()         # get the name of the folder to save the Maxent model results to
        temp = qsdm_settings.getTEMP()+os.sep+"MAXENT" # folder where reprojected files and such are saved
            
        progress.setConsoleInfo("Starting Parameter and File Preperation")
        # Check if temp folder exists, otherwise create it
        if os.path.exists(temp) == False:
            os.mkdir(temp)            
        # Get optional parameters
        param = self.getParameterValue(self.PARAM) 
        o = Processing.getObject(param)        
        if type(o)==QgsVectorLayer and o.isValid() and os.path.splitext(o.source())[1]==".csv":
            progress.setConsoleInfo("Using optional parameter file for MAXENT")

            param = dict()
            # Format Parameters to dictionary
            dp = o.dataProvider()
            for feat in dp.getFeatures():
                geom = feat.geometry()
                com = feat["command"]
                val = feat["value"]
                param[com] = val
            # and make maxent invisible to the modeller
            param["visible"] = False
        else:
            progress.setConsoleInfo("No valid optional Parameter file detected")
            # Use default parameters
            param = dict()
            # per default write a separate maxent results file for each species
            param["perspeciesresults"] = False            
            # and make maxent invisible to the modeller
            param["visible"] = False
        
        # Progress updater:
        n = len(self.getParameterValue(self.ENV).split(";"))+5
        func.updateProcessing(progress,1,n,"Loaded Parameters.")

        
        ## Species layer
        # Get the species file to model. Take selected species column and coordinates
        point = self.getParameterValue(self.SPECIES)
        v = Processing.getObject(point)
        crs = v.crs()
        if crs.authid() != "EPSG:4326":
            progress.setConsoleInfo("Species localities not in WGS84, reprojecting...")
            # Reproject using ogr
            func.reprojectLatLong(v,temp)
            # Then open again as QgsVectorLayer
            out = temp+os.sep+"localities.shp"
            if (os.path.exists(out) and os.path.isfile(out)) == False:
                raise GeoAlgorithmExecutionException("Species point layer data could not be reprojected to WGS84")                 
            fileInfo = QFileInfo(out)
            baseName = fileInfo.baseName()
            v = QgsVectorLayer(out, baseName, "ogr")
            if v.isValid() != True:
                # If this didn't work, try to use the Processing way
                v = Processing.getObject(out)
                if v.isValid() != True:
                    # Otherwise return error
                    raise GeoAlgorithmExecutionException("No valid layer could be loaded from the reprojection.") 
        
        # get names of species from input file
        scl = self.getParameterValue(self.SPEC_COL)                    
        # Get Coordinates from point layer and add the species name
        coord = func.point2table(v,scl)
        if coord is None:
            raise GeoAlgorithmExecutionException("Species point layer data could not be extracted") 

        # Convert coordinates and species name to csv, save in temporary Folder
        # Get Systemwide temporary folder to save the species csv
        speciesPath = temp + os.sep +"species.csv" 
        species = func.saveToCSV(coord,("Species","Long","Lat"),speciesPath)
        specieslist = func.getUniqueAttributeList( v, scl,True)
        progress.setConsoleInfo("Species data successfully prepared for MAXENT")
        progress.setConsoleInfo("---")
        func.updateProcessing(progress,2,n,"Loaded Species data.")
        
        ## Environmental Layers
        # get the selected environmental layers and prepare them for MAXENT 
        progress.setConsoleInfo("Starting preparing the environmental layers")
        envlayers = self.getParameterValue(self.ENV)        
        env = dict()  
        layers = []
        # Project to WGS84 if necessary
        for lay in envlayers.split(";"):
            r = Processing.getObject(lay) # QgsRasterLayer object
            name = str( r.name() )
            crs = r.crs()
            if crs.authid() != "EPSG:4326":
                # Reproject layer
                progress.setConsoleInfo("Originial Layer %s not in WGS84, reprojecting..." % (name))
                r = func.reprojectRasterLatLong(r,temp,True)                
                if r == False or r.isValid()==False :
                    ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Projecting "+name+" to WGS84 failed!")
            layers.append( r.source() )            
        if len(layers) == 0:
            raise GeoAlgorithmExecutionException("Environmental Layers could not be reprojected!") 
        else:
            func.updateProcessing(progress,3,n,"Reprojection finished.")

        # Check the extent of all those layers and unify if necessary
        # Check if necessary -> Do the raster layer have differing extents
        func.updateProcessing(progress,4,n)
        uni = []
        app = False # Which approach should be used?
        if len(layers) > 1 and func.unificationNecessary(layers):
            progress.setConsoleInfo("Input layers have different extents, intersecting...")

            # The credits of the following approach go to Yury Ryabov - http://ssrebelious.blogspot.com
            if app == False:
                # get coordinates of corners for the final raster
                fin_coordinates = func.finCoordinates(layers)
                r = gdal.Open(str( layers[0] ) )
                main_geo_transform = r.GetGeoTransform() 
                proj = r.GetProjection()
                no_data = r.GetRasterBand(1).GetNoDataValue()
                if not no_data:
                    no_data = -9999
                for lay in layers:
                    raster = gdal.Open(str(lay))
                    name = os.path.splitext(os.path.basename(lay))[0]
                    out = temp + os.sep + name + 'warp.tif'                
                    result = func.ExtendRaster(raster, fin_coordinates, out, main_geo_transform, proj, no_data)
                    if result:
                        raster = None
                        if os.path.exists(out):
                            # Add output to uni
                            uni.append(out)
                        else:
                            raise GeoAlgorithmExecutionException("Unified layer could not be saved.") 
                    else:
                        raise GeoAlgorithmExecutionException("Layers could not be unified. Please set do this manually.") 
            else:
                # FIXME: Faster Approach down below. Currently not yet working
                # 1. Build largest extent and geotransform 
                # big_coord has left, top, right, bottom of dataset's bounds in geospatial coordinates.
                fin_coordinates,  main_geo_transform, interp = func.CreateMainGeotransform(layers) # get coordinates and geotransform of corners for the final raster                      
    
                # set number of columns and rows for raster
                main_cols = (fin_coordinates[2] - fin_coordinates[0]) / abs(main_geo_transform[1])  
                main_rows = (fin_coordinates[3] - fin_coordinates[1]) / abs(main_geo_transform[5])
                progress.setConsoleInfo("Creating new raster based on greatest extent with %s columns and %s rows" % (str(main_cols),str(main_rows)))
    
                #FIXME: Check coordinates
                big_coord = [main_geo_transform[0], main_geo_transform[3], main_geo_transform[0] + (main_geo_transform[1] * main_rows), main_geo_transform[3] + (main_geo_transform[5] * main_cols)]
                
                # 2. Loop through rasters and Intersect them export the biggest 
                for lay in layers:
                    name = os.path.splitext(os.path.basename(lay))[0]
                    r = gdal.Open(str( lay ) )
                    src_p = r.GetProjection()
                    if interp:
                        # Interpolate to biggest cellsize
                        progress.setConsoleInfo("Resolution of Environmental Layers is different. Bilinear interpolation to the coarsest cellsize = xy(%s,%s)" % (abs(main_geo_transform[1]),abs(main_geo_transform[5])))
                        #FIXME: Maybe interpolate to nearest neighbor if categorical
                        r = func.gridInterpolation(r,temp,main_geo_transform,main_cols,main_rows,src_p, 'Bilinear',False)
                    wide = abs( r.RasterXSize )
                    high = abs( r.RasterYSize )
                    geotransform = r.GetGeoTransform()
                    nodata = r.GetRasterBand(1).GetNoDataValue() # should be -9999 if projected correctly
                    if nodata == None:
                        nodata = -9999                
                    # target has left, top, right, bottom of dataset's bounds in geospatial coordinates.
                    target = [geotransform[0], geotransform[3], geotransform[0] + (geotransform[1] * wide), geotransform[3] + (geotransform[5] * high)]
                    #Intersection
                    intersection = [max(big_coord[0], target[0]), min(big_coord[1], target[1]), min(big_coord[2], target[2]), max(big_coord[3], target[3])]
                    # Convert to pixels
                    p1 = func.world2Pixel(geotransform,intersection[0],intersection[1])
                    p2 = func.world2Pixel(geotransform,intersection[2],intersection[3])
                    band = r.GetRasterBand(1)            
                    result = band.ReadAsArray(p1[0], p1[1], p2[0] - p1[0], p2[1] - p1[1], p2[0] - p1[0], p2[1] - p1[1])
                    
                    # Write to new raster
                    output = temp + os.sep + name + 'warp.tif'                
                    func.createRaster(output,abs(main_geo_transform[1]),abs(main_geo_transform[5]),result,nodata,main_geo_transform,src_p,'GTiff')
                
                    if os.path.exists(output):
                        # Add output to uni
                        uni.append(output)
                    else:
                        raise GeoAlgorithmExecutionException("Environmental Layers could not be prepared for MAXENT") 
                        return None
        else:
            uni = layers
        if len(uni) == 0 or len(uni) != len(layers):
            raise GeoAlgorithmExecutionException("Environmental Layers with unified extent could not be generated!") 
        else:
            progress.setConsoleInfo("Environmental Layer successfully unified.")
            func.updateProcessing(progress,5,n,"Unified environmental Layers.")
           
        # Format to asc if necessary
        for lay in uni:
            r = Processing.getObject(lay) # QgsRasterLayer object
            name = os.path.basename(str( r.name() ))
            out = temp + os.sep + name + '.asc'
            progress.setConsoleInfo("Convert environmental layers to ESRI ASC format...")
            # Format to asc
            proc = func.raster2ASC(r,out)
            if proc and os.path.isfile(out):
                env[name] = out
            else:
                ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Converting/Projecting "+name+" to ESRI asc format failed!")
        func.updateProcessing(progress,6,n,"Formated to ASC.")
        
        # Check if anything is in env, worked
        if len(env) == 0:
            raise GeoAlgorithmExecutionException("Environmental Layers could not be prepared for MAXENT")
        # Check if the number of the original selected layers is equal to 
        if len(envlayers.split(";")) != len(env):
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"Successfully prepared environmental layers "+str( env.keys() ) )
            raise GeoAlgorithmExecutionException("Not all environmental Layers could be prepared for MAXENT. Check Processing Log.") 
        # Test if species csv exists
        if os.path.exists(speciesPath) == False:
            raise GeoAlgorithmExecutionException("Species point layer could not be prepared for MAXENT") 
        
        ## create the maxent command
        progress.setConsoleInfo("---")
        progress.setConsoleInfo("All fine so far. Attempting to build MAXENT execution command...")

        # Try if JAVA can be executed like this, otherwise take the binary from the given path
        try:
            from subprocess import DEVNULL # python 3k
        except ImportError:
            DEVNULL = open(os.devnull, 'wb')
        proc = subprocess.call(['java', '-version'],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
        if proc == 0:
            start = "java -mx" + str(int(mem)) + "m -jar "
        else:
            progress.setConsoleInfo("JAVA could not be run by default. Using link to binary from set JAVA folder.")
            start = java + " -mx" + str(int(mem)) + "m -jar "

        # if Windows, encapsule jar file in "
        if platform.system() == "Windows":
            start += "\"" + maxent + "\""
        else:
            start += maxent
            
        myCommand = start + " samplesfile=" + speciesPath 
        myCommand += " environmentallayers=" + temp 
                
        # Toggle all selected Layers
        myCommand += " togglelayertype=" 
        for i in range(0,len(env.keys())):
            myCommand += os.path.splitext( env.keys()[i] )[0]
            if i is not len(env.keys())-1:
                myCommand += ","

        myCommand += " outputdirectory=" + work
        # Parse parameters into command
        for option in param.iteritems():
            myCommand += " " + option[0] + "=" + str( option[1] ).lower()
        # finish the command
        myCommand += " redoifexists autorun" 
        
        # add a message
        progress.setConsoleInfo("#### Attempting to start MAXENT ####")
        func.updateProcessing(progress,7,n)

        # execute the command
        loglines = []
        loglines.append('MAXENT execution console output')
#        result = os.system(myCommand)
        proc = subprocess.Popen(
            myCommand,
            shell=True,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            ).stdout
        for line in iter(proc.readline, ''):
            loglines.append(line)
        ProcessingLog.addToLog(ProcessingLog.LOG_INFO, loglines) # Print all loglines if delivered from MAXENT
        err = False
        for line in loglines:
            progress.setConsoleInfo(line)
            if line.find("Error") != -1:
                err = True
        if err:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"MAXENT calculations did not succed! Check the Processing Info output for possible error sources.")
            print "Used command:" + myCommand
        else:
            ProcessingLog.addToLog(ProcessingLog.LOG_INFO,"MAXENT modelling finished.")
            func.updateProcessing(progress,n,n)

            # Finished, Load all resulting layers in QGIS if successfully run
            # In order to be compatible with processing copy or link them to the Processing output folder        
            #out_r = self.getOutputValue(self.OUT_PRED)
            #out_t = self.getOutputValue(self.OUT_PRED_RES)
            p = work + os.sep + "maxentResults.csv"
            func.tableInQgis(p,",")

            #load in only generated Prediction
            for species in specieslist:
                t = species.replace(" ","_")
                p = work + os.sep + t + ".asc"
                func.rasterInQgis(p)
            
            ## Styling and grouping 
            # Freeze the canvas
            canvas = QgsMapCanvas()
            canvas.freeze(True)
            #Add a new group and all new layers to it
            groups = iface.legendInterface().groups() 
            if ('MAXENT' in groups ) == False:
                idx = iface.legendInterface().addGroup( "MAXENT" )
                groups = iface.legendInterface().groups() 
            
            layerMap = QgsMapLayerRegistry.instance().mapLayers()
            for lyr in layerMap.itervalues():                
                if lyr.name() in specieslist:
                    # Move them to the maxent group
                    iface.legendInterface().moveLayer( lyr, groups.index("MAXENT") )                
                    
                    # Style the output
                    lyr.setDrawingStyle("SingleBandPseudoColor")
                    # The band of classLayer
                    classLyrBnd = 1
                    # Color list for ramp
                    clrLst = [  QgsColorRampShader.ColorRampItem(0, QColor(224,224,224),"0"),      # Grey
                                QgsColorRampShader.ColorRampItem(0.01, QColor(0,0,153),"> 0.01"),    # darkblue
                                QgsColorRampShader.ColorRampItem(0.2, QColor(153,204,255),"0.2"),  # lightblue                   
                                QgsColorRampShader.ColorRampItem(0.35,QColor(153,255,153),"0.35"), # lightgreen
                                QgsColorRampShader.ColorRampItem(0.5, QColor(0,153,0),"0.5"),      # green
                                QgsColorRampShader.ColorRampItem(0.65, QColor(255,255,0),"0.65"),  # yellow
                                QgsColorRampShader.ColorRampItem(0.75, QColor(255,128,0),"0.75"),  # orange
                                QgsColorRampShader.ColorRampItem(0.85, QColor(255,0,0),">0.85") ]  # red
                    #Create the shader
                    lyrShdr = QgsRasterShader()
                    #Create the color ramp function
                    clrFnctn = QgsColorRampShader()
                    clrFnctn.setColorRampType(QgsColorRampShader.INTERPOLATED)
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

            #Finally move the Maxent results to the group as well
            lyr = func.getLayerByName( "MaxentResults" )
            iface.legendInterface().moveLayer( lyr, groups.index("MAXENT") )                
            
            canvas.freeze(False)
            canvas.refresh()
    
    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None  
            

class MaxentGUI(GeoAlgorithm):
    
    SPECIES = 'SPECIES'
    SPEC_COL = 'SPEC_COL'

    ENV_DIR = 'ENV_DIR'
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"maxent.png")
    
    def defineCharacteristics(self):
        self.name = 'Maximum Entropy Modelling (Manual Configuration)'
        self.cmdName = 'maxent_gui'
        self.group = 'Species Distribution Modelling'

        self.addParameter(ParameterVector(self.SPECIES, 'Species localities',[ParameterVector.VECTOR_TYPE_POINT,ParameterTable],False)) # Allow point
        self.addParameter(ParameterTableField(Maxent.SPEC_COL, "Species Name", Maxent.SPECIES))
        self.addParameter(ParameterString(self.ENV_DIR, 'Environmental Layers Folder',qsdm_settings.workPath() ))

    
    def processAlgorithm(self, progress):
        # Get the location to the maxent jar file
        maxent = qsdm_settings.maxent()
        if os.path.basename(maxent) == 'maxent.jar':
            pass
        else:
            maxent = os.path.join(qsdm_settings.maxent(), 'maxent.jar') # If the directory and not the file was chosen
        # Get location of java, available memory and output path
        if sys.platform == "win32" or "win64":
            ex = "java.exe"
        else:
            ex = "java"
        java = os.path.join(qsdm_settings.javaPath(),ex) # the path to java if basic execution fails
        temp = qsdm_settings.getTEMP()+os.sep+"MAXENT" # folder where reprojected files and such are saved
        mem = str( qsdm_settings.getMEM() )     # available memory for MAXENT
        work = qsdm_settings.workPath()         # get the name of the folder to save the Maxent model results to        
        env_dir = self.getParameterValue(self.ENV_DIR)
        progress.setConsoleInfo("Starting Species Layer Preperation")
        # Check if temp folder exists, otherwise create it
        if os.path.exists(temp) == False:
            os.mkdir(temp)
        ## Species layer preperation
        # Get the species file to model. Take selected species column and coordinates
        point = self.getParameterValue(self.SPECIES)
        v = Processing.getObject(point)
        scl = self.getParameterValue(self.SPEC_COL) # get names of species from input file
        if v.source().find("type=csv") != -1 :
            raise GeoAlgorithmExecutionException("Species point layer should be saved as ESRI Shapefile")                 
        else:
            crs = v.crs()
            if crs.authid() != "EPSG:4326":
                progress.setConsoleInfo("Species localities not in WGS84, reprojecting...")
                # Reproject using ogr
                func.reprojectLatLong(v,temp)
                # Then open again as QgsVectorLayer
                out = temp+os.sep+"localities.shp"
                if (os.path.exists(out) and os.path.isfile(out)) == False:
                    raise GeoAlgorithmExecutionException("Species point layer data could not be reprojected to WGS84")                 
                fileInfo = QFileInfo(out)
                baseName = fileInfo.baseName()
                v = QgsVectorLayer(out, baseName, "ogr")
                if v.isValid() != True:
                    # If this didn't work, try to use the Processing way
                    v = Processing.getObject(out)
                    if v.isValid() != True:
                        # Otherwise return error
                        raise GeoAlgorithmExecutionException("No valid layer could be loaded from the reprojection.") 
            
            # Get Coordinates from point layer and add the species name
            coord = func.point2table(v,scl)
            if coord is None:
                raise GeoAlgorithmExecutionException("Species point layer data could not be extracted") 
    
            # Convert coordinates and species name to csv, save in temporary Folder
            # Get Systemwide temporary folder to save the species csv
            speciesPath = temp + os.sep +"species.csv" 
            species = func.saveToCSV(coord,("Species","Long","Lat"),speciesPath)
            specieslist = func.getUniqueAttributeList( v, scl, True)
        progress.setConsoleInfo("Species data successfully prepared for MAXENT")
        progress.setConsoleInfo("---")
        ## Maxent execution
        # Try if JAVA can be executed like this, otherwise take the binary from the given path
        try:
            from subprocess import DEVNULL # python 3k
        except ImportError:
            DEVNULL = open(os.devnull, 'wb')
        proc = subprocess.call(['java', '-version'],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
        if proc == 0:
            start = "java -mx" + str(int(mem)) + "m -jar "
        else:
            progress.setConsoleInfo("JAVA could not be run by default. Using link to binary from set JAVA folder.")
            start = java + " -mx" + str(int(mem)) + "m -jar "

        # if Windows, encapsule jar file in "
        if platform.system() == "Windows":
            start += "\"" + maxent + "\""
        else:
            start += maxent
            
        myCommand = start + " samplesfile=" + speciesPath 
        myCommand += " environmentallayers=" + env_dir 

        myCommand += " outputdirectory=" + work
        # finish the command
        myCommand += " redoifexists" 
        
        # add a message
        progress.setConsoleInfo("#### Attempting to start MAXENT ####")

        # execute the command
        loglines = []
        loglines.append('MAXENT execution console output')
#        result = os.system(myCommand)
        proc = subprocess.Popen(
            myCommand,
            shell=True,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            ).stdout
        for line in iter(proc.readline, ''):
            loglines.append(line)
        ProcessingLog.addToLog(ProcessingLog.LOG_INFO, loglines) # Print all loglines if delivered from MAXENT
        err = False
        for line in loglines:
            progress.setConsoleInfo(line)
            if line.find("Error") != -1:
                err = True
        if err:
            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,"MAXENT calculations did not succed! Check the Processing Info output for possible error sources.")
            print "Used command:" + myCommand
        else:
            ProcessingLog.addToLog(ProcessingLog.LOG_INFO,"MAXENT modelling finished.")

            # Finished, Load all resulting layers in QGIS if successfully run
            # In order to be compatible with processing copy or link them to the Processing output folder        
            #out_r = self.getOutputValue(self.OUT_PRED)
            #out_t = self.getOutputValue(self.OUT_PRED_RES)
            p = work + os.sep + "maxentResults.csv"
            func.tableInQgis(p,",")

            #load in only generated Prediction
            for species in specieslist:
                t = species.replace(" ","_")
                p = work + os.sep + t + ".asc"
                func.rasterInQgis(p)
                # Finished

    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None

class MaxentParameters(GeoAlgorithm):
    
    OUT_F = 'OUT_F'
    of = ['logistic','raw','cumulative']
    RESP = 'RESP'
    RESP_EXP = 'RESP_EXP'
    PIC = 'PIC'
    JACK = 'JACK'
    RANDOM = 'RANDOM'
    LOGSCALE = 'LOGSCALE'
    CLAMPGRID = 'CLAMPGRID'
    MESS = 'MESS'
    RANDOM_POINTS = 'RANDOM_POINTS'
    BETA_MULT = 'BETA_MULT'
    REPLICATES = 'REPLICATES'
    REPLICATE_TYPE = 'REPLICATE_TYPE'
    R_TYPE = ["Crossvalidate","Bootstrap","Subsample"]
    LINEAR = 'LINEAR'
    QUADRATIC = 'QUADRATIC'
    PRODUCT = 'PRODUCT'
    THRESHOLD = 'THRESHOLD'
    HINGE = 'HINGE'
    FADEBYCLAMPING = 'FADEBYCLAMPING'
    EXTRAPOLATE = 'EXTRAPOLATE'
    PLOTS = 'PLOTS'
    MAXITERATIONS = 'MAXITERATIONS'
    CONVG_THRESH = 'CONVG_THRESH'
    PROC= 'PROC'
    DEF_PREV = 'DEF_PREV'
    APPLY_THRESH = 'APPLY_THRESH'
    PERSPECRES = 'PERSPECRES'
    
    # Out file
    OUT_PARAM = 'OUT_PARAM'
    
    def getIcon(self):
        return QIcon(os.path.dirname(__file__) +os.sep+".."+ os.sep+"icons"+os.sep+"maxent.png")
    
    def defineCharacteristics(self):
        self.name = 'MAXENT (ParameterPreperation)'
        self.cmdName = 'maxent_param'
        self.group = 'Species Distribution Modelling'

        self.addParameter(ParameterSelection(self.OUT_F, "Output Format", self.of, 0))
        self.addParameter(ParameterBoolean(self.RESP, 'Create Response Curves',False))
        self.addParameter(ParameterBoolean(self.RESP_EXP, 'Show exponent in Response Curves',False))
        self.addParameter(ParameterBoolean(self.PIC, 'Create pictures for each grid',True))
        self.addParameter(ParameterBoolean(self.JACK, 'Measure Importance by jackknifing',False))
        self.addParameter(ParameterBoolean(self.RANDOM, 'Random seed',False))
        self.addParameter(ParameterBoolean(self.LOGSCALE, 'Display Pictures on a logscale',False))
        self.addParameter(ParameterBoolean(self.CLAMPGRID, 'Write Grid of spatial Clamping',True))
        self.addParameter(ParameterBoolean(self.MESS, 'Write MESS grid',True))
        self.addParameter(ParameterNumber(self.RANDOM_POINTS, 'Percentage of presence localities to be randomly set aside as test points',0,None,0))
        self.addParameter(ParameterNumber(self.BETA_MULT, 'Multiply all automatic regularization parameters by this number.',1.0,None,1.0))
        self.addParameter(ParameterNumber(self.REPLICATES, 'Number of Replicates for crossvalidation',False,True,1))
        self.addParameter(ParameterSelection(self.REPLICATE_TYPE, "What Type of Replicate", self.R_TYPE, 0))
        self.addParameter(ParameterBoolean(self.LINEAR, 'Allow linear features to be used',True))
        self.addParameter(ParameterBoolean(self.QUADRATIC, 'Allow quadratic features to be used',True))
        self.addParameter(ParameterBoolean(self.PRODUCT, 'Allow product features to be used',True))
        self.addParameter(ParameterBoolean(self.THRESHOLD, 'Allow threshold features to be used',True))
        self.addParameter(ParameterBoolean(self.HINGE, 'Allow hinge features to be used',True))
        self.addParameter(ParameterBoolean(self.FADEBYCLAMPING, 'Reduce prediction by the difference between clamped and non-clamped output',False))
        self.addParameter(ParameterBoolean(self.EXTRAPOLATE, 'Predicts to regions out of environmental space',True))
        self.addParameter(ParameterBoolean(self.PLOTS, 'Write various plots for inclusion in .html output',True))
        self.addParameter(ParameterNumber(self.MAXITERATIONS, 'Stop training after this many iterations of the optimization algorithm',1,None,500))
        self.addParameter(ParameterNumber(self.CONVG_THRESH, 'Stop training when the drop in log loss per iteration drops below this number',False,True,0.00001))
        self.addParameter(ParameterNumber(self.PROC, 'Number of processor threads to use',False,True,1))
        self.addParameter(ParameterNumber(self.DEF_PREV, 'Default prevalence of the species: probability of presence at ordinary occurrence points.',False,True,0.5))
#        self.addParameter(ParameterSelection(self.APPLY_THRESH, "Apply a threshold rule",""))
        self.addParameter(ParameterNumber(self.CONVG_THRESH, 'Stop training when the drop in log loss per iteration drops below this number',False,True,0.00001))
        self.addParameter(ParameterBoolean(self.PERSPECRES, 'Create Per Species results',True))  
              
        self.addOutput(OutputTable(self.OUT_PARAM,'MAXENT_Parameters'))                

    def processAlgorithm(self, progress):
        param = dict()
        param["outputformat"] = self.of[self.getParameterValue(self.OUT_F)]
        param["responsecurves"] = self.getParameterValue(self.RESP)
        param["responsecurvesexponent"] = self.getParameterValue(self.RESP_EXP)
        param["pictures"] = self.getParameterValue(self.PIC)
        param["randomseed"] = self.getParameterValue(self.RANDOM)
        param["logscale"] = self.getParameterValue(self.LOGSCALE)
        param["writeclampgrid"] = self.getParameterValue(self.CLAMPGRID)
        param["writemess"] = self.getParameterValue(self.MESS)
        param["randomtestpoints"] = self.getParameterValue(self.RANDOM_POINTS)
        param["betamultiplier"] = self.getParameterValue(self.BETA_MULT)
        param["replicates"] = self.getParameterValue(self.REPLICATES)
        param["replicatetype"] = self.R_TYPE[self.getParameterValue(self.REPLICATE_TYPE)]
        param["linear"] = self.getParameterValue(self.LINEAR)
        param["quadratic"] = self.getParameterValue(self.QUADRATIC)
        param["product"] = self.getParameterValue(self.PRODUCT)
        param["threshold"] = self.getParameterValue(self.THRESHOLD)
        param["hinge"] = self.getParameterValue(self.HINGE)
        param["fadebyclamping"] = self.getParameterValue(self.FADEBYCLAMPING)
        param["extrapolate"] = self.getParameterValue(self.EXTRAPOLATE)
        param["plots"] = self.getParameterValue(self.PLOTS)
        param["maximumiterations"] = self.getParameterValue(self.MAXITERATIONS)
        param["convergencethreshold"] = self.getParameterValue(self.CONVG_THRESH)
        param["threads"] = self.getParameterValue(self.PROC)
        param["defaultprevalence"] = self.getParameterValue(self.DEF_PREV)
        param["perspeciesresults"] = self.getParameterValue(self.PERSPECRES)           
        
#        param["applythresholdrule"]= self.getParameterValue(self.APPLY_THRESH)
        
        res = []
        for option in param.iteritems():
            res.append((option[0],str( option[1] ).lower()))
        
        out = self.getOutputValue(self.OUT_PARAM)
        species = func.saveToCSV(res,("command","value"),out)

    def help(self):
        helppath = os.path.join(os.path.dirname(__file__) + os.sep + ".." + os.sep + "help", self.cmdName + ".html")
        if os.path.isfile(helppath):
            return False, helppath
        else:
            return False, None