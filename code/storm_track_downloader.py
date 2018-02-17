# -*- coding: utf-8 -*-
import urllib2,ssl
ssl._create_default_https_context = ssl._create_unverified_context
import urllib
import base64
import os
from bs4 import BeautifulSoup #For extracting links from nhc html
import re #For refining extracted links to those we actually want
import requests
import getpass
import logging

#This function downloads storm forecasts from multiple sources over HTML
#PARAMETERS:
#     runFilepath - the rundir of the simulation these storm forecasts will belong to
#     fileSource - specify which forecast source you wish to use 
def DownloadStorms(runFilepath, fileSource, targetYear, username, password):
    dir = runFilepath + '/input/storm_forecasts/'
    FORECAST_DIR = dir
    fileList =[]
    downloadCounter = 0

    if fileSource == 'nhc':
        logging.info("*** Beginning NHC Forecast Download ***")
        REMOTE_BASE_URL = 'https://manati.star.nesdis.noaa.gov/TC_cone_info/'

        #Get directory page for forecasts
        try:
            response = urllib2.urlopen(REMOTE_BASE_URL)
        except (urllib2.HTTPError, urllib2.URLError), err:
            StormErrorHandler(fileSource, err)
            return #Abort attempted storm download

        soup = BeautifulSoup(response, "lxml")

        for link in soup.findAll('a', attrs={'href': re.compile("cone_info.txt$")}):
            fileList.append(link.get('href'))

        #Download all files (to ensure newest version)
        for file in fileList:
            #Find year associated with forecast
            date = file.split('_')[0]
            year = date[4:8]

            url = REMOTE_BASE_URL + file
            filename = FORECAST_DIR + file

            if year == targetYear:
                #print file
                try:
                    urllib.urlretrieve(url, filename)
                    downloadCounter += 1
                except Exception, err:
                    StormErrorHandler(fileSource, err)
                    return #Abort further storm downloads

        logging.info("%s NHC storm(s) downloaded.", downloadCounter)
        downloadCounter = 0

    elif fileSource == 'jtwc':
        logging.info("*** Beginning JTWC Forecast Download ***")
        REMOTE_BASE_URL = 'https://pzal.ndbc.noaa.gov/atcf_storms/storms/'
        
        #Must set up opener since site is password protected
        #username = raw_input('Username for https://pzal.ndbc.noaa.gov: ');
        #password = getpass.getpass()
        p = urllib2.HTTPPasswordMgrWithDefaultRealm()
        p.add_password(None, REMOTE_BASE_URL, username, password)
        handler = urllib2.HTTPBasicAuthHandler(p)
        opener = urllib2.build_opener(handler)
        urllib2.install_opener(opener)

        #Get directory page for forecasts
        try:
            response = urllib2.urlopen(REMOTE_BASE_URL)
        except (urllib2.HTTPError, urllib2.URLError), err:
            StormErrorHandler(fileSource, err)
            return #Abort attempted storm download

        soup = BeautifulSoup(response, "lxml")
    
        for link in soup.findAll('a', attrs={'href': re.compile(".fst$")}):
            fileList.append(link.get('href'))

        for file in fileList:
            url = REMOTE_BASE_URL + file
            filename = FORECAST_DIR + file

            year = file[4:8]

            if year == targetYear:
                try:
                    r = requests.get(url, auth=(username, password))
                    write_file = open(filename, "w")
                    write_file.write(r.text)
                    write_file.close()

                    downloadCounter +=1
                except Exception, err:
                    StormErrorHandler(fileSource, err)
                    return #Abort further storm downloads

        logging.info("%s JTWC storm(s) downloaded.", downloadCounter)

    else:
        logging.error("ERROR: Invalid Storm Forecast Source %s. Skipping download", fileSource)
        return

def StormErrorHandler(source, err):
    if type(err).__name__  == 'urllib2.HTTPError':
        logging.error("HTTPError %s %d", err.reason, err.code)
    else:
        logging.error("URLError %s", err)

    logging.error("%s forecast download aborted - some or all files may not have been downloaded", source)



    
