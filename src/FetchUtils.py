# AA2CODON - A Python package for generating codon MSA file from amino acid MSA file
  #   Copyright (C) 2015  Etai Jacob
  #   etai.jacob@gmail.com
  #
  #   This program is free software; you can redistribute it and/or modify
  #   it under the terms of the GNU General Public License as published by
  #   the Free Software Foundation; either version 2 of the License, or
  #   (at your option) any later version.
  #
  #   This program is distributed in the hope that it will be useful,
  #   but WITHOUT ANY WARRANTY; without even the implied warranty of
  #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #   GNU General Public License for more details.
  #
  #   You should have received a copy of the GNU General Public License along
  #   with this program; if not, write to the Free Software Foundation, Inc.,
  #   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



# WSDL URL for service.
wsdlUrl = 'http://www.ebi.ac.uk/ws/services/WSDbfetchDoclit?wsdl'

# Libs for downloading utilities:
import platform
import os
import sys
import urllib2
import logging

import suds
from suds.client import Client

__author__ = 'etai'

# Setup logging
logging.basicConfig(level=logging.INFO)

# Output level
outputLevel = 1
# Debug level
debugLevel = 0

def getChunks(l, n):
    return [l[i:i + n] for i in range(0, len(l), n)]


# Fetch a set of entries.
def soapFetchBatch(dbName, idListStr, formatName, styleName):
    entriesStr = dbfetch.fetchBatch(dbName, idListStr, formatName, styleName)
    return entriesStr


# Print an entry.
def printFetchData(query, formatName, styleName):
    entryStr = soapFetchData(query, formatName, styleName)
    print entryStr


# Print a set of entries.
def returnFetchBatch(dbName, idListStr, formatName, styleName):
    entriesStr = soapFetchBatch(dbName, idListStr, formatName, styleName)
    return entriesStr


# Debug print
def printDebugMessage(functionName, message, level):
    if (level <= debugLevel):
        print >> sys.stderr, '[' + functionName + '] ' + message


def soap_setup():
    # Create the service interface
    printDebugMessage('main', 'WSDL: ' + wsdlUrl, 1)

    client = Client(wsdlUrl)
    if outputLevel > 1:
        print client
    global dbfetch
    dbfetch = client.service

    # Set the client user-agent.
    clientRevision = '$Revision: 2467 $'
    clientVersion = '0'
    if len(clientRevision) > 11:
        clientVersion = clientRevision[11:-2]
    userAgent = 'EBI-Sample-Client/%s (%s; Python %s; %s) suds/%s Python-urllib/%s' % (
        clientVersion, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        suds.__version__, urllib2.__version__
    )
    printDebugMessage('main', 'userAgent: ' + userAgent, 1)
    httpHeaders = {'User-agent': userAgent}
    client.set_options(headers=httpHeaders)

    # Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
    proxyOpts = dict()
    if os.environ.has_key('http_proxy'):
        proxyOpts['http'] = os.environ['http_proxy'].replace('http://', '')
    elif os.environ.has_key('HTTP_PROXY'):
        proxyOpts['http'] = os.environ['HTTP_PROXY'].replace('http://', '')
    if 'http' in proxyOpts:
        client.set_options(proxy=proxyOpts)
