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
import numpy
# Try to import functions from osgeo
try:
    from osgeo import gdal, gdalconst
except ImportError:
    import gdal
try:
    from osgeo import ogr, osr
except ImportError:
    import ogr,osr

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.ProcessingLog import ProcessingLog
from processing.core.Processing import Processing
from processing.core.ProcessingConfig import ProcessingConfig
from processing.core.GeoAlgorithmExecutionException import GeoAlgorithmExecutionException
    
# Converts a raster to ASC using gdal
# Expects a QgsRasterLayer returned by Processing
def raster2ASC(raster,out):
    rasterPath = str( raster.source() )
    # Get nodata
    srcImage = gdal.Open(str( rasterPath ))
    band = srcImage.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if nodata == None:
        #set nodata to -9999 per default
        nodata = -9999
    nodata = str(nodata)
    # Use the gdal_translate command
    #FIXME: Maybe look for a within python solution
    try:
        from subprocess import DEVNULL # python 3k
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')
    proc = subprocess.call(['gdal_translate', '-of', 'AAIGrid','-a_nodata',nodata,rasterPath, out],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
    if proc == 0:        
        return True
    else:
        return False
    
# Reprojects a raster using gdalwarp
# Expects a QgsRasterLayer from Processing and returns one if successful as well
def reprojectRasterLatLong(raster,targetDir,console=True):
    if console == False:
        
        rasterPath = str( raster.source() )
        srcImage = gdal.Open(rasterPath)
        driver = srcImage.GetDriver()
        band = srcImage.GetRasterBand(1)
        data_type = band.DataType

        gt= srcImage.GetGeoTransform()
        cols = srcImage.RasterXSize
        rows = srcImage.RasterYSize
        ext = GetExtentFromGT(gt,cols,rows)
        
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt( srcImage.GetProjection() )
        tgt_srs = osr.SpatialReference()
        tgt_srs.ImportFromEPSG(4326)
        tgt_srs = src_srs.CloneGeogCS()
        
        geo_ext = ReprojectCoords(ext,src_srs,tgt_srs) # Coorner Coordinates in lat-long
        cols = (geo_ext[2] - geo_ext[0]) / gt[1]  
        rows = (geo_ext[3] - geo_ext[1]) / gt[5]

        metadata = driver.GetMetadata()
        if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
            pass
        else:
            progress.setConsoleInfo("Creation of input Fileformat is not supported by gdal. Create GTiff by default.")
            driver = gdal.GetDriverByName("GTiff")            
        
        out = targetDir + os.sep + name + '_WGS84.tif'
        try:
            outData = driver.Create(out, cols, rows, 1, data_type)
        except Exception, e:
            raise GeoAlgorithmExecutionException("Output reprojected Raster file could not be created!")
        
        band = outData.GetRasterBand(1)
        band.WriteArray( band.ReadAsArray() )
        band.FlushCache()
        na = srcImage.GetRasterBand(1).GetNoDataValue()
        if na is not None:
            band.SetNoDataValue( na )
        else:
            band.SetNoDataValue( -9999 )

        new_gt = [geo_ext[0],cols,gt[2],geo_ext[3],gt[4],rows]
        outData.SetGeoTransform(new_gt) 
        outData.SetProjection(tgt_srs)
        outData = None # Close writing

        fileInfo = QFileInfo(out)
        baseName = fileInfo.baseName()
        rlayer = QgsRasterLayer(out, baseName)
        return rlayer
        
    else:
        rasterPath = str( raster.source() )
        # Get nodata
        srcImage = gdal.Open(str( rasterPath ))
        #FIXME: Easier way in gdal to get Projection EPSG?
        s_crs = raster.dataProvider().crs()
        proj = s_crs.authid()
        band = srcImage.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        if nodata == None:
            #set nodata to -9999 per default
            nodata = -9999
        nodata = str(nodata)
        name = str( raster.name() )
        out = targetDir + os.sep + name + '_WGS84.tif'
        try:
            from subprocess import DEVNULL # python 3k
        except ImportError:
            DEVNULL = open(os.devnull, 'wb')
        proc = subprocess.call(['gdalwarp','-overwrite', '-of', 'GTiff','-s_srs',proj,'-t_srs','EPSG:4326',rasterPath, out],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
        if proc == 0:
            fileInfo = QFileInfo(out)
            baseName = fileInfo.baseName()
            rlayer = QgsRasterLayer(out, baseName)
            return rlayer
        else:
            return False

def GetExtentFromGT(gt,cols,rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]

    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
        yarr.reverse()
    return ext

def InvGeoTransform(gt_in):
  '''
   ************************************************************************
   *                        InvGeoTransform(gt_in)
   ************************************************************************

   **
   * Invert Geotransform.
   *
   * This function will invert a standard 3x2 set of GeoTransform coefficients.
   *
   * @param  gt_in  Input geotransform (six doubles - unaltered).
   * @return gt_out Output geotransform (six doubles - updated) on success,
   *                None if the equation is uninvertable.
  '''
  #    ******************************************************************************
  #    * This code ported from GDALInvGeoTransform() in gdaltransformer.cpp
  #    * as it isn't exposed in the python SWIG bindings until GDAL 1.7
  #    * copyright & permission notices included below as per conditions.
  #
  #    ******************************************************************************
  #    * $Id: gdaltransformer.cpp 15024 2008-07-24 19:25:06Z rouault $
  #    *
  #    * Project:  Mapinfo Image Warper
  #    * Purpose:  Implementation of one or more GDALTrasformerFunc types, including
  #    *           the GenImgProj (general image reprojector) transformer.
  #    * Author:   Frank Warmerdam, warmerdam@pobox.com
  #    *
  #    ******************************************************************************
  #    * Copyright (c) 2002, i3 - information integration and imaging
  #    *                          Fort Collin, CO
  #    *
  #    * Permission is hereby granted, free of charge, to any person obtaining a
  #    * copy of this software and associated documentation files (the "Software"),
  #    * to deal in the Software without restriction, including without limitation
  #    * the rights to use, copy, modify, merge, publish, distribute, sublicense,
  #    * and/or sell copies of the Software, and to permit persons to whom the
  #    * Software is furnished to do so, subject to the following conditions:
  #    *
  #    * The above copyright notice and this permission notice shall be included
  #    * in all copies or substantial portions of the Software.
  #    *
  #    * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
  #    * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  #    * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  #    * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  #    * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  #    * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  #    * DEALINGS IN THE SOFTWARE.
  #    ****************************************************************************

  # we assume a 3rd row that is [1 0 0]

  # Compute determinate
  det = gt_in[1] * gt_in[5] - gt_in[2] * gt_in[4]

  if( abs(det) < 0.000000000000001 ):
      return

  inv_det = 1.0 / det

  # compute adjoint, and divide by determinate
  gt_out = [0,0,0,0,0,0]
  gt_out[1] =  gt_in[5] * inv_det
  gt_out[4] = -gt_in[4] * inv_det

  gt_out[2] = -gt_in[2] * inv_det
  gt_out[5] =  gt_in[1] * inv_det

  gt_out[0] = ( gt_in[2] * gt_in[3] - gt_in[0] * gt_in[5]) * inv_det
  gt_out[3] = (-gt_in[1] * gt_in[3] + gt_in[0] * gt_in[4]) * inv_det

  return gt_out


def ReprojectCoords(coords,src_srs,tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates
    '''
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords

# Reproject shapefile to WGS84 using ogr
def reprojectLatLong(layer,target):
    targetRef = osr.SpatialReference()
    targetRef.ImportFromEPSG(4326)         # from EPSG
    
    vectorPath = layer.source()

    # Create an OGR layer from a boundary shapefile used to clip
    shapef = ogr.Open("%s" % str(vectorPath))
    lyr = shapef.GetLayer()
    sourceRef = lyr.GetSpatialRef()
        
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(sourceRef, targetRef)
    
    # create the output layer
    driver = ogr.GetDriverByName('ESRI Shapefile')
    outputShapefile = target + os.sep + "localities.shp"
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)
    outLayer = outDataSet.CreateLayer("species_4326", geom_type=ogr.wkbPoint)

    # add fields
    inLayerDefn = lyr.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = lyr.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = lyr.GetNextFeature()
    
    # Create the prj file
#     spr = outLayer.GetSpatialRef()
#     spr.MorphToESRI()
#     file = open(target + 'localities.prj', 'w')
#     file.write(spatialRef.ExportToWkt())
#     file.close()
    
    # close the shapefiles
    shapef.Destroy()
    outDataSet.Destroy()

# General function to retrieve layers
def getLayerByName( layerName ):
  layerMap = QgsMapLayerRegistry.instance().mapLayers()
  for name, layer in layerMap.iteritems():
    if layer.name() == layerName:
        if layer.isValid():
          return layer
        else:
          return None
          
# Get all field values of a given attribute from a vector layer
def getUniqueAttributeList( vlayer, field):
  path = vlayer.source()
  datasource = ogr.Open(str(path))
  layer = datasource.GetLayer()
  layerName = layer.GetName()
  field = str(field)
  try:
    d = datasource.ExecuteSQL("SELECT %s FROM %s" % (field,layerName))
  except RuntimeError:
    QMessageBox.warning(QDialog(),"LecoS: Warning","Failed to query the vector layers attribute table")
    return
  attr = []
  for i in range(0,d.GetFeatureCount()):
    f = d.GetFeature(i)
    val = f.GetField(0)
    val = val.replace(" ","_") # Replace spaces with underscores
    if val not in attr:
        attr.append(val)
  return attr
        
# Exracts point coordinates to a table
def point2table(layer,scl):
    coord = []
    if type(layer) == QgsVectorLayer and layer.isValid():
        dp = layer.dataProvider()
        for feat in dp.getFeatures():
            geom = feat.geometry()
            name = feat[str( scl )]
            coord.append(( name,geom.asPoint().x(),geom.asPoint().y() )) 
        return coord
    else:
        return None


# Saves table to csv
def saveToCSV(results, titles, filePath ):
    f = open(filePath, "wb" )
    writer = csv.writer(f,delimiter=',',quotechar="",quoting=csv.QUOTE_NONE)
    writer.writerow(titles)
    for item in results:
        writer.writerow(item)
    f.close()
    
# Adds a generated Raster to the QGis table of contents
def rasterInQgis(rasterPath):
  fileName = str(rasterPath)
  fileInfo = QFileInfo(fileName)
  baseName = fileInfo.baseName()
  rlayer = QgsRasterLayer(fileName, baseName)
  if not rlayer.isValid():
    raise GeoAlgorithmExecutionException("Layer is not valid. Failed to add the generated Layer to QGis") 
  QgsMapLayerRegistry.instance().addMapLayer(rlayer)
  

# Adds a vector layer to the QGis table of contents
def tableInQgis(vectorPath,delim):
  fileName = str(vectorPath)
  fileInfo = QFileInfo(fileName)
  baseName = fileInfo.baseName()
  uri = "file:/"+fileName+"?delimiter=%s" % (delim)
  vlayer = QgsVectorLayer(uri, baseName, "delimitedtext")
  if not vlayer.isValid():
    raise GeoAlgorithmExecutionException("LecoS: Warning","Failed to add the Layer to QGis")
  QgsMapLayerRegistry.instance().addMapLayer(vlayer)

# Create basic raster without projection
def createRaster(output,cols,rows,array,nodata,gt,proj=None,d='GTiff'):
    driver = gdal.GetDriverByName(d)
    # Create File based in path
    try:
        tDs = driver.Create(output, cols, rows, 1, gdal.GDT_Float32)
    except Exception, RuntimeError:
        raise GeoAlgorithmExecutionException("Could not generate output file")
    
    # set the NoData value
    band = tDs.GetRasterBand(1)
    band.WriteArray(array)

    # georeference the image and set the projection
    tDs.SetGeoTransform(gt)

    if proj != None:
        try:
            tDs.SetProjection(proj)
        except Exception, err:
            pass

    try:
      band.SetNoDataValue(nodata)
    except TypeError:
      pass#band.SetNoDataValue(-9999) # set -9999 in the meantime
    
    # flush data to disk
    band.FlushCache()
    band = tDs = None # Close writing

def getPixelSize(rasterPath):
    """
    Takes a path as input
    """
    fileName = str(rasterPath)
    fileInfo = QFileInfo(fileName)
    baseName = fileInfo.baseName()
    rlayer = QgsRasterLayer(fileName, baseName)
    x = rlayer.rasterUnitsPerPixelX() # Extract The X-Value
    y = rlayer.rasterUnitsPerPixelY() # Extract The X-Value
    return x,y

def getArrayFromRaster(rasterPath):
    try:
        raster = gdal.Open(rasterPath)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()
        return array        
    except Exception, AttributeError:
        return None

def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate 
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line)    

def Pixel2world(self,geoMatrix, x, y):
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    coorX = (ulX + (x * xDist))
    coorY = (ulY + (y * yDist))
    return (coorX, coorY)

def ApplyGeoTransform(inx,iny,gt):
  ''' Apply a geotransform
      @param  inx:       Input x coordinate (double)
      @param  iny:       Input y coordinate (double)
      @param  gt:        Input geotransform (six doubles)

      @return: outx,outy Output coordinates (two doubles)
  '''
  # Copyright (c) 2013 Australian Government, Department of the Environment
  #
  # Permission is hereby granted, free of charge, to any person obtaining a copy
  # of this software and associated documentation files (the "Software"), to deal
  # in the Software without restriction, including without limitation the rights
  # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  # copies of the Software, and to permit persons to whom the Software is
  # furnished to do so, subject to the following conditions:
  #
  # The above copyright notice and this permission notice shall be included in
  # all copies or substantial portions of the Software.
  #
  # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  # THE SOFTWARE.

  outx = gt[0] + inx*gt[1] + iny*gt[2]
  outy = gt[3] + inx*gt[4] + iny*gt[5]
  return (outx,outy)

def MapToPixel(mx,my,gt):
  ''' Convert map to pixel coordinates
      @param  mx:    Input map x coordinate (double)
      @param  my:    Input map y coordinate (double)
      @param  gt:    Input geotransform (six doubles)
      @return: px,py Output coordinates (two ints)

      @change: changed int(p[x,y]+0.5) to int(p[x,y]) as per http://lists.osgeo.org/pipermail/gdal-dev/2010-June/024956.html
      @change: return floats
      @note:   0,0 is UL corner of UL pixel, 0.5,0.5 is centre of UL pixel
  '''
  # Copyright (c) 2013 Australian Government, Department of the Environment
  #
  # Permission is hereby granted, free of charge, to any person obtaining a copy
  # of this software and associated documentation files (the "Software"), to deal
  # in the Software without restriction, including without limitation the rights
  # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  # copies of the Software, and to permit persons to whom the Software is
  # furnished to do so, subject to the following conditions:
  #
  # The above copyright notice and this permission notice shall be included in
  # all copies or substantial portions of the Software.
  #
  # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  # THE SOFTWARE.

  if gt[2]+gt[4]==0: #Simple calc, no inversion required
      px = (mx - gt[0]) / gt[1]
      py = (my - gt[3]) / gt[5]
  else:
      px,py=ApplyGeoTransform(mx,my,InvGeoTransform(gt))
  #return int(px),int(py)
  return px,py

def PixelToMap(px,py,gt):
  ''' Convert pixel to map coordinates
      @param  px:    Input pixel x coordinate (double)
      @param  py:    Input pixel y coordinate (double)
      @param  gt:    Input geotransform (six doubles)
      @return: mx,my Output coordinates (two doubles)

      @note:   0,0 is UL corner of UL pixel, 0.5,0.5 is centre of UL pixel
  '''
  # Copyright (c) 2013 Australian Government, Department of the Environment
  #
  # Permission is hereby granted, free of charge, to any person obtaining a copy
  # of this software and associated documentation files (the "Software"), to deal
  # in the Software without restriction, including without limitation the rights
  # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  # copies of the Software, and to permit persons to whom the Software is
  # furnished to do so, subject to the following conditions:
  #
  # The above copyright notice and this permission notice shall be included in
  # all copies or substantial portions of the Software.
  #
  # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  # THE SOFTWA
  mx,my=ApplyGeoTransform(px,py,gt)
  return mx,my

def ExtendRaster(raster, xy_list, output, main_geo_transform, proj, no_data):
  '''
  Extends canvas of the given raster to the extent provided by xy_list in
  [minx, miny, maxx, maxy] format

  @param raster   Input raster to be processed (GDALDataset)
  @param xy_list  Input list of min and max XY coordinates of resulting raster
                  (list)
  @param output   Input full path to the resulting raster
  @param main_geo_transform: Main parameters of geotransformation (gdal.GetGeoTransform())
  @param proj:    Projection of the raster (gdal.GetProjection())
  @param no_data: a value to be set as 'No Data' value
  @returm         Output indicates that function ended successfully (only True)
  '''
  # check xy_list
  if len(xy_list) != 4:
    sys.exit('xy_list must contain 4 values!')

  minx = xy_list[0]
  miny = xy_list[1]
  maxx = xy_list[2]
  maxy = xy_list[3]

  band = raster.GetRasterBand(1)
  NA = band.GetNoDataValue()
  data_type = band.DataType

  # set number of columns and rows for raster
  geo_transform = raster.GetGeoTransform()
  columns = (maxx - minx) / main_geo_transform[1]
  columns = int(abs(columns))
  rows = (maxy - miny) / main_geo_transform[5]
  rows = int(abs(rows))
  bands = raster.RasterCount

  # create raster
  r_format = "GTiff"
  driver = gdal.GetDriverByName(r_format)
  metadata = driver.GetMetadata()
  if metadata.has_key( gdal.DCAP_CREATE ) and metadata[ gdal.DCAP_CREATE ] == "YES":
    pass
  else:
    print "Driver %s does not support Create() method." % format
    return False

  try:
    outData = driver.Create(output, columns, rows, bands, data_type)
  except:
    return NoPathGiven()
  outData.SetProjection(proj)

  # we don't want rotated raster in output
  new_geo_transform = [minx, main_geo_transform[1], 0.0,
    maxy, 0.0, main_geo_transform[5]]
  outData.SetGeoTransform(new_geo_transform)

  for i in xrange(1, (bands +1) ):
    band = raster.GetRasterBand(i)
    band = band.ReadAsArray()

    # we create array here to write into it created raster
    new_raster = numpy.zeros( (rows, columns) )

    # populate array with the values from original raster
    for col in xrange(columns):
      for row in xrange(rows):
        # convert array value location to XY-coordinates
        x,y = PixelToMap(col, row, new_geo_transform)

        # convert coordinates to pixel location at the original raster
        px,py = MapToPixel(x, y, geo_transform)

        # extract pixel value
        try:
          if px < 0 or py < 0:
            pix_value = no_data
          else:
            pix_value = band[py, px]
            if pix_value == NA:
              pix_value = no_data
        except:
          pix_value = no_data

        # assign extracted value to array
        new_raster[row, col] = pix_value

    outData.GetRasterBand(i).WriteArray(new_raster)

  # close dataset
  outData = None

  return True



def gridInterpolation(raster,temp,match_geotrans,main_cols,main_rows, match_proj=None, interp='Bilinear',command = False):
    '''
    Expects a gdal raster layer as input
    '''
    if command == False:
        if raster == None:
            raise GeoAlgorithmExecutionException("Raster could not be interpolated. Something before went wrong.")         
        # Processing
        name = os.path.splitext(os.path.basename( raster.GetDescription()) )[0]
        output = temp + os.sep + name + ".tif"
        src_proj = raster.GetProjection()
        src_geotrans = raster.GetGeoTransform()
        nodata = raster.GetRasterBand(1).GetNoDataValue() # keep the nodata value
        if nodata == None:
            nodata = -9999
        wide = raster.RasterXSize
        high = raster.RasterYSize 
        # Output / destination
        try:
            # try create File driver based in path
            dst = gdal.GetDriverByName('GTiff').Create(output, wide, high, 1, gdalconst.GDT_Float32)
        except RuntimeError:
            raise GeoAlgorithmExecutionException("Could not generate interpolated output file")
            
        dst.SetGeoTransform( match_geotrans )
        dst.SetProjection( match_proj)
        dst.GetRasterBand(1).SetNoDataValue(int( nodata ) ) # write old nodata value
    
        # Do the work
        if interp == 'Bilinear':
            gdal.ReprojectImage(raster, dst, src_proj, match_proj, gdalconst.GRA_Bilinear)
        elif interp == 'Cubic':
            gdal.ReprojectImage(raster, dst, src_proj, match_proj, gdalconst.GRA_Cubic)        
        elif interp == 'Cubicspline':
            gdal.ReprojectImage(raster, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)            
        elif interp == 'Lanczos':
            gdal.ReprojectImage(raster, dst, src_proj, match_proj, gdalconst.GRA_Lanczos)            
        elif interp == 'NearestNeighbour':
            gdal.ReprojectImage(raster, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
                    
        return dst
    else:
        # Using the gdalwarp command-line for interpolation
        # gdalwarp -te xmin ymin xmax ymax -tr 3 3 -r bilinear A_state_30m.tif C_county_3m.tif
        # where:
        #     -te target extents, you need to supply this from your (B) county-sized raster; try using gdalinfo to help determine this extent
        #     -tr target resolution, 3 m
        #     -r bilinear, a good algorithm for orthophotos and DEMs, but not for others; other algorithms are available
        #     A_state_30m.tif, input (A) file
        #     C_county_3m.tif, output (C) file
        rasterPath = str( raster.GetDescription() )
        name = os.path.splitext(os.path.basename( raster.GetDescription()) )[0]
        output = str( temp + os.sep + name + "warp.tif" )
        nodata = raster.GetRasterBand(1).GetNoDataValue() # keep the nodata value
        if nodata == None:
            nodata = -9999 #otherwise set to -9999
        nodata = str(nodata)
        ex = FindCorners(raster)
        # Keep resolution
        x,y = getPixelSize(rasterPath)    
        # Try and see if file already exists
        if os.path.exists(output):
            # Remove warped temporary file
            try:
                os.remove(output)
            except Exception, e:
                raise GeoAlgorithmExecutionException("Error creating warped raster in temporary folder.")
        try:
            from subprocess import DEVNULL # python 3k
        except ImportError:
            DEVNULL = open(os.devnull, 'wb')
        proc = subprocess.call(['gdalwarp','-te',str( ex[0] ),str( ex[1] ),str( ex[2] ),str( ex[3] ), '-tr',str(x),str(y),'-r',interp.lower(), rasterPath, output],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
        print ['gdalwarp','-te',str( ex[0] ),str( ex[1] ),str( ex[2] ),str( ex[3] ), '-tr',str(x),str(y),'-r',interp.lower(), rasterPath, output]
        if proc == 0:
            r = gdal.Open(output)
            if r is None:
                raise GeoAlgorithmExecutionException("Error calling gdalwarp to interpolate to maximal extent.")
            else:
                return r
        else:
            raise GeoAlgorithmExecutionException("Error calling gdalwarp for interpolation.")

#def vectorToRaster(self, fieldName, layerName, layer, referenceRasterName, outputRaster):
# 		# band set
# 		if self.bndSetPresent == "Yes" and referenceRasterName == self.bndSetNm:
# 			referenceRasterName = self.bndSet[0]
# 			# input
# 			refRstr = self.selectLayerbyName(referenceRasterName)
# 		else:
# 			if self.selectLayerbyName(referenceRasterName) is None:
# 				self.msg4()
# 				self.refreshRasterLayer()
# 			else:
# 				# input
# 				refRstr = self.selectLayerbyName(referenceRasterName)
# 		# register drivers
# 		gdal.AllRegister()
# 		try:
# 			refRstrSrc = refRstr.source().encode(self.fSEnc)
# 			rstrCheck = "Yes"
# 		except Exception, err:
# 			# logger
# 			if self.logSttngVal == "Yes": self.logToFile(str(inspect.stack()[0][3])+ " " + self.lineOfCode(), " ERROR exception: " + str(err))
# 			rstrCheck = "No"
# 		if rstrCheck == "No":
# 			self.msg4()
# 		else:
# 			# open input with GDAL
# 			refRstrDt = gdal.Open(refRstrSrc, GA_ReadOnly)
# 			# number of x pixels
# 			refRstrCols = refRstrDt.RasterXSize
# 			# number of y pixels
# 			refRstrRows = refRstrDt.RasterYSize
# 			# check projections
# 			refRstrProj = refRstrDt.GetProjection()
# 			# pixel size and origin
# 			refRstGeoTrnsf = refRstrDt.GetGeoTransform()
# 			refRstPxlXSz = abs(refRstGeoTrnsf[1])
# 			refRstPxlYSz = abs(refRstGeoTrnsf[5])
# 			tifDrvr = gdal.GetDriverByName( "GTiff" )
# 			outputRaster = tifDrvr.Create(outputRaster, refRstrCols, refRstrRows, 1, GDT_Int32)
# 			outputRasterBand = outputRaster.GetRasterBand(1)
# 			# set raster projection from reference
# 			outputRaster.SetGeoTransform( [ refRstGeoTrnsf[0] , refRstGeoTrnsf[1] , 0 , refRstGeoTrnsf[3] , 0 , refRstGeoTrnsf[5] ] )
# 			outputRaster.SetProjection(refRstrProj)
# 			outputRasterBand.SetNoDataValue(-9999)
# 			matrix = numpy.zeros((refRstrRows, refRstrCols), dtype='int32')
# 			matrix.fill(-9999)
# 			outputRasterBand.WriteArray(matrix, 0, 0)
# 			outputRasterBand.FlushCache()
# 			source_ds = ogr.Open(layer)
# 			source_layer = source_ds.GetLayer()
# 			# convert reference layer to raster
# 			outCheck = gdal.RasterizeLayer(outputRaster, [1], source_layer, options = ["ATTRIBUTE=" + str(fieldName)])
# 			# close bands
# 			outputRasterBand = None
# 			# close rasters
# 			outputRaster = None
# 			# logger
# 			if self.logSttngVal == "Yes": self.logToFile(str(inspect.stack()[0][3])+ " " + self.lineOfCode(), "vector to raster check: " + outCheck)


def unificationNecessary(rasterList):
    ex = None
    for ras in rasterList:
        r = gdal.Open(ras)
        ck = FindCorners(r)
        if ex == None:
            ex = ck
        else:
            if ck != ex:
                return True
        r = None
    return False

def FindCorners(raster):
    '''
    returns min and max X and Y values of the given raster
    
    @param raster:               Input raster (GDALDataset)
    @return list of coordinates  Output min and max coordinates of raster (list)
    '''
    width = raster.RasterXSize
    height = raster.RasterYSize
    geo_transform = raster.GetGeoTransform()
    top_left_x = geo_transform[0]
    top_left_y = geo_transform[3]
    top_right_x = geo_transform[0] + width*geo_transform[1]
    top_right_y = geo_transform[3] + width*geo_transform[4]
    bottom_right_y = geo_transform[3] + width*geo_transform[4] + height*geo_transform[5]
    bottom_right_x = geo_transform[0] + width*geo_transform[1] + height*geo_transform[2]
    bottom_left_x =  geo_transform[0] + 1*geo_transform[1] + height*geo_transform[2]
    bottom_left_y = geo_transform[3] + height*geo_transform[5]
    
    x_list = [top_left_x, top_right_x, bottom_right_x, bottom_left_x]
    y_list = [top_left_y, top_right_y, bottom_left_y, bottom_right_y]
    
    min_x = min(x_list)
    max_x = max(x_list)
    min_y = min(y_list)
    max_y = max(y_list)
    
    return [min_x, min_y, max_x, max_y]
        
def finCoordinates(rasters_list):
    '''
    returns a list of the coordinates for the unified raster
    
    @param rasters_list: Lists of rasters (list)
    @return:             List of coordinates [min_x, min_y, max_x, max_y]
    '''
    fin_coordinates = []
    
    for raster in rasters_list:
        rast = gdal.Open(raster)
        coordinates = FindCorners(rast)
        rast = None
        if len(fin_coordinates) != 4:
            fin_coordinates = coordinates
        else:
            minx = coordinates[0]
            miny = coordinates[1]
            maxx = coordinates[2]
            maxy = coordinates[3]
            if minx < fin_coordinates[0]:
                fin_coordinates[0] = minx
            if miny < fin_coordinates[1]:
                fin_coordinates[1] = miny
            if maxx > fin_coordinates[2]:
                fin_coordinates[2] = maxx
            if maxy > fin_coordinates[3]:
                fin_coordinates[3] = maxy
    
    return fin_coordinates

def CreateMainGeotransform(rasters_list):
    """
    gt[0] /* top left x */
    gt[1] /* w-e pixel resolution */ --> Pixelwidth
    gt[2] /* rotation, 0 if image is "north up" */
    gt[3] /* top left y */
    gt[4] /* rotation, 0 if image is "north up" */
    gt[5] /* n-s pixel resolution */ --> Pixelheight
    """
    fin_coordinates = finCoordinates(rasters_list)
    # Try to see if Pixelresolution is equal and test if rotation is equal everywhere
    width = []
    height = []
    ori = []
    for raster in rasters_list:
        src = gdal.Open(raster)
        width.append( src.RasterXSize )
        height.append( src.RasterYSize )
        ori.append( (src.GetGeoTransform()[2],src.GetGeoTransform()[4] ) )
        src = None
    if len( numpy.unique(width) ) > 1 or len( numpy.unique(height) ) > 1:
        ProcessingLog.addToLog(ProcessingLog.LOG_WARNING,"QSDM: Environmental layers cellsize were not equal. Interpolating to the biggest common cellsize. Consider doing this beforehand.")
        interp = True
        width = max(width) # inverse max
        height = min(height)*-1 # inverse min
    else:
        width = width[0]
        height = height[0]
        interp = False
    if numpy.sum(ori) != 0:
        ProcessingLog.addToLog(ProcessingLog.LOG_WARNING,"QSDM: Some of the rasterlayers don't have a north up orientation. Results may be flawed")
    else:
        ori = [0.0,0.0] # Northwards
        
    gt = [fin_coordinates[0],width,ori[0],fin_coordinates[3],ori[1],height]
    return fin_coordinates, gt, interp
    
# Update Processing Process
def updateProcessing(progress,i,n,text=None):
    progress.setPercentage(int(100 * i / n))
    if text != None:
        progress.setText(text)
    