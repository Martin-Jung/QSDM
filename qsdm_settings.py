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

import os,sys, tempfile
from qgis.core import QgsApplication
import subprocess

from processing.core.ProcessingConfig import ProcessingConfig
from processing.core.ProcessingLog import ProcessingLog
from processing.tools.system import *


class qsdm_settings:

    JAVA_EXEC = 'JAVA_EXEC'
    MAXENT = 'MAXENT'
    WORK_DIR = 'WORK_DIR'
    MEM = 'MEM'
    TEMP = 'TEMP'

    @staticmethod
    def javaPath():
        # load previous value
        folder = ProcessingConfig.getSetting(qsdm_settings.JAVA_EXEC)
        if folder is None:
            # Try to automatically detect JAVA path
            # Can Java be executed from the shell?            
            try:
                from subprocess import DEVNULL # python 3k
            except ImportError:
                DEVNULL = open(os.devnull, 'wb')
            try:
                proc = subprocess.call(['java', '-version'],stdin=subprocess.PIPE, stdout=DEVNULL, stderr=subprocess.STDOUT)
            except Exception, err:
                ProcessingLog.addToLog(ProcessingLog.LOG_ERROR, "There has been an Error detecting JAVA for QSDM. Please set correct the path manually.")   
                folder = ''
                proc = 1 
            if proc == 0:
                # java detected, look for executable
                env = os.environ
                #a = env.get("PATH") bin also in here
                if env.has_key("JAVA_HOME"):
                    folder = os.path.join(env.get("JAVA_HOME"), 'bin')
            else:
                # Otherwise User has to set path manually
                folder = ''
        return folder

    @staticmethod
    def maxent():
        folder = ProcessingConfig.getSetting(qsdm_settings.MAXENT)
        if folder is None:
            folder = ''
        return folder

    @staticmethod
    def workPath():
        folder = ProcessingConfig.getSetting(qsdm_settings.WORK_DIR)
        if folder is None:
            folder = tempfile.gettempdir()
        return folder
    
    @staticmethod
    def getMEM():
        mem = ProcessingConfig.getSetting(qsdm_settings.MEM)
        return mem
        
    @staticmethod
    def getTEMP():
        temp = ProcessingConfig.getSetting(qsdm_settings.TEMP)
        if temp is None:
            temp = tempfile.gettempdir()
        return temp

