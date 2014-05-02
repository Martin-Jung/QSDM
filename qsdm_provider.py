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
from PyQt4.QtGui import *

# Processing bindings
from processing.core.AlgorithmProvider import AlgorithmProvider
from processing.core.ProcessingConfig import ProcessingConfig
from processing.core.ProcessingConfig import Setting
from processing.core.ProcessingLog import ProcessingLog

# Custom Processing Settings file
from qsdm_settings import qsdm_settings

# Import modules
import os

# The main algorithm Provider. Modules will be added to Processing (SEXTANTE) 
# if their requirements are met

class qsdm_serviceProvider(AlgorithmProvider):  
    def __init__(self):
        AlgorithmProvider.__init__(self)   
        # Check requirements
        self.checkRequirements()        
        # Create algorithms list
        self.createAlgsList()
    
    def getDescription(self):
        return "QSDM (species distribution modelling)"

    def initializeSettings(self):
        '''In this method we add settings needed to configure our provider.
        Do not forget to call the parent method, since it takes care or
        automatically adding a setting for activating or deactivating the
        algorithms in the provider'''
        AlgorithmProvider.initializeSettings(self)
        
        # Path to Java Jar
        ProcessingConfig.addSetting(Setting(self.getDescription(),
                                    qsdm_settings.JAVA_EXEC,
                                    'Detected JAVA folder with java executable',
                                    qsdm_settings.javaPath()))
        # Path to MAXENT
        ProcessingConfig.addSetting(Setting(self.getDescription(),
                                    qsdm_settings.MAXENT,
                                    'Path to maxent.jar file',
                                    qsdm_settings.maxent()))

        # Working folder
        ProcessingConfig.addSetting(Setting(self.getDescription(),
                                    qsdm_settings.WORK_DIR,
                                    'MAXENT Working Folder',
                                    qsdm_settings.workPath()))
        
        # Memory for SDM
        ProcessingConfig.addSetting(Setting(self.getDescription(),
                                    qsdm_settings.MEM,
                                    'Available Memory for Computation',
                                    '512'))
        # Temporary folder
        ProcessingConfig.addSetting(Setting(self.getDescription(),
                                    qsdm_settings.TEMP,
                                    'Temporary Folder',
                                    qsdm_settings.getTEMP()))
        
        # R path

        
        
        '''To get the parameter of a setting parameter, use SextanteConfig.getSetting(name_of_parameter)'''

    def unload(self):
        '''Setting should be removed here, so they do not appear anymore
        when the plugin is unloaded'''
        AlgorithmProvider.unload(self)
        ProcessingConfig.removeSetting(qsdm_settings.JAVA_EXEC)
        ProcessingConfig.removeSetting(qsdm_settings.MAXENT)
        ProcessingConfig.removeSetting(qsdm_settings.MEM)        
        ProcessingConfig.removeSetting(qsdm_settings.WORK_DIR)
        ProcessingConfig.removeSetting(qsdm_settings.TEMP)

    def getName(self):
        '''This is the name that will appear on the toolbox group.
        It is also used to create the command line name of all the algorithms
        from this provider'''
        return "qsdm"

    def getIcon(self):
        '''We return the icon for qsdm'''
        return QIcon(os.path.dirname(__file__) + os.sep + "icons"+os.sep+"default.png")

    def checkRequirements(self):
        self.req = {}
        # Check for JAVA
        if qsdm_settings.javaPath() == '' or None:
            self.req["JAVA"] = False
        else:
            self.req["JAVA"] = True
        # Check for MAXENT
        if qsdm_settings.maxent() == '' or None or len(qsdm_settings.maxent()) < 4:
            self.req["MAXENT"] = False
        else:
            self.req["MAXENT"] = True        
        # Check for Numpy
        try:
            import numpy as np
            self.req["Numpy"] = True
        except ImportError:
            self.req["Numpy"] = False        
        # Check for Scipy
        try:
            import scipy as sp
            self.req["Scipy"] = True
        except ImportError:
            self.req["Scipy"] = False
        # Check for Scikits
        self.req["Scikits"] = False
        # R support available?
        try:
            import rpy2
            self.req["R"] = True
        except ImportError:
            self.req["R"] = False            
        
#         # ---- Create Log outputs ---- #
#         ProcessingLog.addToLog(ProcessingLog.LOG_INFO,
#                     "QSDM:JAVA found: " + str(self.req["JAVA"]))
# 
#         ProcessingLog.addToLog(ProcessingLog.LOG_INFO,
#                     "QSDM:Scipy found: " + str(self.req["Scipy"]))
            
    def createAlgsList(self):
        '''Create list of Arguments based on system-wide configuration and available libraries'''
        
        self.preloadedAlgs = []        
        # Check available libraries and load tools accordingly
        
        # Data Preperation methods require only numpy
        if self.req["Numpy"] == True:
            from algorithms.DataPreperation import *
            #self.preloadedAlgs.append( VectorOutlierSelection() )
            self.preloadedAlgs.append( CreateRichnessGrid() )
            self.preloadedAlgs.append( DataTransformationSimple() )
        
        # Maxent support if JAVA is running and enable the tool. Give Error message to inform user to 
        # specify the path to the maxent binary if not already set
        if self.req["JAVA"] == True:
            from algorithms.Maxent import *
            self.preloadedAlgs.append( MaxentParameters() ) # Generate optional table displaying the MAXENT parameters
            self.preloadedAlgs.append( Maxent() ) # Add Maxent Modelling
            self.preloadedAlgs.append( MaxentGUI() ) # Add Maxent Modelling
         
#         if self.req["Scipy"] == True:
#             try:
#                 import algorithms.LinearRegression

        for alg in self.preloadedAlgs:
            alg.provider = self # reset provider

    
    def _loadAlgorithms(self):
        '''Here we fill the list of algorithms in self.algs.
        This method is called whenever the list of algorithms should be updated.
        If the list of algorithms can change while executing SEXTANTE for QGIS
        (for instance, if it contains algorithms from user-defined scripts and
        a new script might have been added), you should create the list again
        here.
        In this case, since the list is always the same, we assign from the pre-made list.
        This assignment has to be done in this method even if the list does not change,
        since the self.algs list is cleared before calling this method'''
        self.algs = self.preloadedAlgs