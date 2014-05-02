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
 This script initializes the plugin, making it known to QGIS.
"""
def name():
    return "QSDM"
def description():
    return "Species distribution modelling support for the QGIS Processing toolbox"
def version():
    return "Version 0.1"
def icon():
    return "icons/default.png"
def qgisMinimumVersion():
    return "2.0"
def classFactory(iface):
    from qsdm_main import qsdm
    return qsdm(iface)
