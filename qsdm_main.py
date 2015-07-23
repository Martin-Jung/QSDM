# -*- coding: utf-8 -*-
"""
/***************************************************************************
 QSDM
<<<<<<< HEAD
        Species distribution modelling support for the QGIS Processing toolbox
=======
        Species distribution modelling support for the QGIS Processing toolb    ox
>>>>>>> old
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

from PyQt4 import QtGui
# Import QGIS analysis tools in case we need them
from qgis.core import *
from qgis.gui import *
from qgis.analysis import *

# Import base libraries
import os,re,sys,csv,string,math,operator,subprocess,tempfile,inspect

# Try to import Sextante
try:
    from processing.core.Processing import Processing
    sex_load = True
except ImportError:
    sex_load = False

if sex_load:
    # Add folder to sys.path
    cmd_folder = os.path.split(inspect.getfile( inspect.currentframe() ))[0]
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)

# Import QSDM provider
from qsdm_provider import qsdm_serviceProvider

## Basic Processing Class ##
class qsdm( object ):
    def __init__(self, iface):
        # Save reference to the QGIS interface
        self.iface = iface
                
        # Link to Provider
        self.provider = qsdm_serviceProvider()        

    # Init Provider
    def initGui(self):
        Processing.addProvider(self.provider)

    # Unload provider
    def unload(self):
        Processing.removeProvider(self.provider)