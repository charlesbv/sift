#replace this line with the path to Canopy python"


# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import matplotlib as mpl
mpl.use('WXAgg')
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.animation as animation
from collections import *
import time
import datetime
from datetime import timedelta
from itertools import groupby
from operator import itemgetter
from math import *
from scipy.interpolate import *
from bisect import *
import urllib
import matplotlib.gridspec as gridspec
import logging
import ctypes
import wx.lib.scrolledpanel
#Non-Library function imports
from formatMapString import *
from read_input_file import *
from find_in_read_input_order_variables import *
from storm_track_downloader import DownloadStorms
from readZoom import readZoom
from spock_local import *
import cygIcon #For importing window/taskbar icon


import glob
from pprint import pprint
import os
import pandas as pd
from scipy.spatial import distance
import numpy as np

# for cbv: #!/Users/cbv/Library/Enthought/Canopy_64bit/User/bin/python

#****************************************************

try:
    import wx
    import wx.lib.colourdb as cdb #For color keys
except ImportError:
    raise ImportError, "The wxPython module is required to run this application."

# if not sys.platform.startswith('win'):
#     pathToApp = __file__
#     while not pathToApp.endswith('code'):
#         pathToApp = os.path.dirname( pathToApp )
#     os.chdir(pathToApp)


#Set up logging - level=logging.INFO is default, for additional debug messages change to level=logging.DEBUG
logging.basicConfig(filename='../sift.log', level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout)) #Makes logging also print to stdout
sys.stderr=open('../sift.log', 'a+') #This allows all crashes to be dumped to log file

#Print log header
logging.info("*** Starting SIFT ***")



demo_names = ['short', 'long', 'storm']
#Non-standard running modes
if len(sys.argv) > 1:
    #Demo Mode
    if sys.argv[1] == "demo" or sys.argv[1] == "Demo":
        if len(sys.argv) != 3:
            logging.error('Demo argument error - incorrect number of arguments')
            print "Demo mode requires an additional argument specifying which demo to run.\n" \
            + "Please enter a command of the form: ./sift.py demo <demo_name> where\n" \
            + "demo name is one of:"
            print demo_names
            sys.exit()

        if sys.argv[2] == 'short':
            startDate = '2017-04-02T00:00:00'
            endDate = '2017-04-02T06:00:00'
        elif sys.argv[2] == 'long':
            startDate = '2017-03-22T00:00:00'
            endDate = '2017-03-24T00:00:00'
        elif sys.argv[2] == 'storm':
            startDate = '2017-02-10T00:00:00'
            endDate = '2017-02-10T12:00:00'
        else:
            logging.error("Invalid demo name - %s", sys.argv[2])
            print "Please enter a command of the form: ./sift.py demo <demo_name> where\n" \
            + "demo name is one of:"
            print demo_names
            sys.exit()

        runFilepath = (startDate + "_end_" + endDate).replace(":", "_")
        inputFilepath = '../demo_input/' + runFilepath

    #Help command
    elif sys.argv[1] == "help" or sys.argv[1] == "Help":
        print "For help, please refer to the SIFT user manual available in the sift/doc directory."
        sys.exit()

#Default running mode
else: 
    # Get analysis start and end date
    # v1.2
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    startDate = raw_input('Analysis Start Date (YYYY-MM-DDTHH:MM:SS): ')
    endDate = raw_input('Analysis Stop Date (YYYY-MM-DDTHH:MM:SS): ')
#     startDate = '2018-02-01T00:00:00'
#     endDate = '2018-02-01T00:05:00'
    if ( 'T' in startDate ) == False:
        startDate = startDate + "T00:00:00"
    if ( 'T' in endDate ) == False:
        endDate = endDate + "T00:00:00"

    # Set filepath for rundir of this simulation
    runFilepath = (startDate + "_end_" + endDate).replace(":", "_")
    inputFilepath = '../input_sift/' + runFilepath
    
    # ask for password to download storm files 
    # v1.2
    username_jtwc = raw_input('Username for https://pzal.ndbc.noaa.gov: '); 
    if sys.stdin.isatty():
        password_jtwc = getpass.getpass('Password: ')
    else:
        print "Password (it'll show on the screen): "
        password_jtwc =  sys.stdin.readline().rstrip() 


    # # Run SpOCK propagator to generate necessary SIFT inputs
    
    if not os.path.isdir(inputFilepath):
        #logging.info("***Starting CYGNSS Server Connection***")
        spock_local(startDate, endDate)

    # print "DONE"
    # raise Exception
    #Assumes the year you are targeting is the year you are looking at
    targetYear = startDate[0:4]
    
    # Download storm forecasts for this simulation
    # v1.2
    DownloadStorms(inputFilepath, 'nhc', targetYear, username_jtwc, password_jtwc)
    DownloadStorms(inputFilepath, 'jtwc', targetYear, username_jtwc, password_jtwc) #Note: requires login, to suppress login request simply comment out

logging.info("*** Starting main SIFT execution (this may take a few minutes) ***")
driverFile = "spock_spec_start_" + startDate.replace(":","_") + "_end_" + endDate.replace(":","_")
logging.debug("%s",inputFilepath)


class SatelliteManager(object):
    animDirection = 0 #Integer that indicates forward vs backward animation; 1 = forward, -1 = backward, 0 = paused
    curSample = 0
    speedFactor = 1 #Integer that determines speed of animation
    firstTime = ''
    lastTime = ''



    satColors = ['lawngreen', 'blue', 'purple', 'mediumorchid', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
    colorByIntensity = {'DB':'deepskyblue', 'TD':'darkturquoise', 'TS':'orange', 'TY':'darkcyan', 'ST':'darkslategray', 'TC':'tomato',
                                    'HU':'red', 'SD':'chartreuse', 'SS':'forestgreen', 'EX':'purple', 'PT':'limegreen',
                                    'IN':'blueviolet', 'DS':'darkorchid', 'LO': 'slateblue', 'WV': 'mediumslateblue', 'ET':'blueviolet',
                                    'XX':'slategray'}
    #This dictionary is meant to hold a 'lighter' version of the main colors representing each storm type for forecasts in future
    colorByIntensity_faded = {'DB':'skyblue', 'TD':'paleturquoise', 'TS':'navajowhite', 'TY':'lightcyan', 'ST':'slategray', 'TC':'lightsalmon',
                                    'HU':'pink', 'SD':'palegreen', 'SS':'lightgreen', 'EX':'orchid', 'PT':'yellowgreen',
                                    'IN':'thistle', 'DS':'orchid', 'LO': 'azure', 'WV': 'lightblue', 'ET':'plum',
                                    'XX':'whitesmoke'}



    #User dictated file locations - ASSUMPTION: sift.py located in ./code
    driver_filepath = inputFilepath + '/input/spock_in/' + driverFile +'.txt'
    storm_dir = inputFilepath + '/input/storm_forecasts'
    sat_filepath = inputFilepath + '/input/sat_positions/interpolated_position_LLA_ECEF_'
    specular_filepath = inputFilepath + '/input/spec_positions/specular_'
    specInteraction_filepath = specular_filepath
    gsInteraction_filepath = inputFilepath + '/input/gsCoverage/report_all_by_'
    screenshot_filepath = '../output_sift/' + runFilepath +'/screenshots'
    gsReportLoc = '../output_sift/' + runFilepath + '/cygnss_gs_reports'
    specReportLoc = '../output_sift/' + runFilepath +'/cygnss_spec_reports'
    overpassReportLoc = '../output_sift/' + runFilepath + '/cygnss_overpass_reports'


    def __init__(self, m, maxlonIn, minlonIn, latIn, heightIn, widthIn, axMapIn, inZoom):

        self.step_visualization = 1
        self.nb_storms = 0
        self.stormsShowed = False
        self.inZoom = inZoom

        #These parameters are used only for placing clock on map, copied from map instantiation
        self.max_lon = maxlonIn
        self.min_lon = minlonIn
        self.min_lat = latIn
        self.height_fig_with_map = heightIn
        self.width_fig_with_map = widthIn
        self.ax_map = axMapIn

        #IMPORTANT NOTE - SATELLITE NUMBERING
        #Due to the nature of the CYGNSS TLEs
        # - sc 1 is FM05                                                                                                       
        # - sc 2 is FM04                                                                                                       
        # - sc 3 is FM02                                                                                                       
        # - sc 4 is FM01                                                                                                      
        # - sc 5 is FM08                                                                                                       
        # - sc 6 is FM06                                                                                                       
        # - sc 7 is FM07                                                                                                       
        # - sc 8 is FM03

        ###Accumulate satellite data###
        self.satNames = ['CYGFM05', 'CYGFM04', 'CYGFM02', 'CYGFM01', 'CYGFM08', 'CYGFM06', 'CYGFM07', 'CYGFM03']
        self.specNames = ['SPEC1', 'SPEC2', 'SPEC3', 'SPEC4']

        #Read propagator input files
        self.input_variables, self.order_input_variables = read_input_file(self.driver_filepath)
        self.date_start = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'date_start')]
        self.date_stop = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'date_stop')]
        self.dt = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'dt')]
        self.nb_steps = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'nb_steps')]
        self.nb_satellites = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'nb_sc')]
        self.gps_name = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'gps_name')]
        self.output_file_path_propagator = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'output_file_path_list')]
        self.output_filename_propagator = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'output_file_name_list')]
        #self.nb_storms = self.input_variables[9]
        self.storm_name = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'storm_name')]

        #Read propagator output files
        self.nb_spec_pts = 4
        self.interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
        self.nb_steps_interpolation = (int)((self.nb_steps-1) * self.dt / self.interpolation_step) +1

        #Declare necessary variables
        self.lon_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
        self.lat_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
        self.ecef_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation, 3])
        self.heading_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
        self.time_sat = []
        self.lon_spec, self.lat_spec, self.gain_spec, self.distance_spec_to_storm, self.specular_over_storm = np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation])
        self.ecef_spec = np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation, 3])
        self.name_spec = []
        self.time_spec = []
        self.date = []
        self.list_output_variables_to_read = ["longitude","latitude"]
        self.point = namedtuple('point', ['x', 'y'])
        self.color = namedtuple('color', 'red green blue')
        self.spacecraft_list = []
        self.specular_list = []
        self.init_index = 3 #27, 5, WHO KNOWS

        #Build out satellite data

        #SAT TRACE LINES - used to show ground tracks
        self.satTraceX = [[] for i in range(self.nb_satellites)]
        self.satTraceY = [[] for i in range(self.nb_satellites)]
        self.specTraceX = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceY = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.satTraceRefs = {}
        self.specTraceRefs = {}
        self.specGainTraceRefs = {}

        #Ground station interaction parameters
        self.GSinteractionIdx = None

        for i in range(self.nb_satellites):
            self.time_spec_sublist = []
            self.name_spec_between_list_and_sublist = []

            ### SATELLITES ###
            output_file = open(self.sat_filepath + self.output_filename_propagator[i], "r")
            read_output_file = output_file.readlines()
            nb_lines_header_output_file_sat = 0
            
            for j in range(self.nb_steps_interpolation):

                    
                
                try:
                    self.lon_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[1])
                    if (self.lon_sat[i,j] > 180):
                        self.lon_sat[i,j] = self.lon_sat[i,j] - 360.
                    self.lat_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[2])
                    self.ecef_sat[i,j,0] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[4])
                    self.ecef_sat[i,j,1] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[5])
                    self.ecef_sat[i,j,2] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[6])
                    self.heading_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[7])
                    self.time_sat.append(read_output_file[j+nb_lines_header_output_file_sat].split()[0])
                except:
                    # print read_output_file[j+nb_lines_header_output_file_sat]
                    # print 'XXX'
                    # print self.sat_filepath + self.output_filename_propagator[i]
                    # raise Exception
                    if j > len(read_output_file)-nb_lines_header_output_file_sat-1:#this can happen if SpOCK stopped outputing before the end (sometimes it does that one time step before the end)
                        self.lon_sat[i,j] = self.lon_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1]
                        if (self.lon_sat[i,j] > 180):
                            self.lon_sat[i,j] = self.lon_sat[i,j] - 360.
                        self.lat_sat[i,j] = self.lat_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1]
                        self.ecef_sat[i,j,0] = self.ecef_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1,0]
                        self.ecef_sat[i,j,1] = self.ecef_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1,1]
                        self.ecef_sat[i,j,2] = self.ecef_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1,2]
                        self.heading_sat[i,j] = self.heading_sat[i,len(read_output_file)-nb_lines_header_output_file_sat-1]
                        time_temp_str = self.time_sat[-1]
                        time_temp = datetime.strptime(time_temp_str, "%Y-%m-%dT%H:%M:%S")
                        time_new = time_temp + timedelta(seconds = self.interpolation_step )
                        time_new_str = datetime.strftime(time_new, "%Y-%m-%dT%H:%M:%S")
                        self.time_sat.append(time_new_str)
                        #print self.time_sat[-1]
                    else:                        
                        print 'Trouble with this one.'

            # Build the tuples for the visualization of the satellites
            spacecraft = namedtuple('spacecraft',('name',) +  self.point._fields + ('point_plot',) + ('marker_spacecraft',))
            self.spacecraft_list.append(spacecraft)
            name_temp = self.output_filename_propagator[i].replace(".txt","")
            self.spacecraft_list[i].name = name_temp
            # initial position
            self.spacecraft_list[i].x, self.spacecraft_list[i].y =  m(self.lon_sat[i,self.init_index], self.lat_sat[i,self.init_index])
            self.spacecraft_list[i].marker_spacecraft = 'v'
            # point on the plot
            self.spacecraft_list[i].point_plot = m.plot([],[],  marker=self.spacecraft_list[i].marker_spacecraft, markersize=10,color = self.satColors[i])[0]

            ### SPECULAR POINTS ###
            self.spec_dir = ""
            for j in range(len(self.output_file_path_propagator[i].split('/'))-2):
                if (j > 0):
                    self.spec_dir = self.spec_dir + "/" + self.output_file_path_propagator[i].split('/')[j]

            file_specular = open(self.specular_filepath + self.output_filename_propagator[i], "r")
            read_file_specular  = file_specular.readlines()
            # Nb of lines in the spec file header
            if (i == 0):
                nb_lines_header_output_file_spec = 0
                while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
                    nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
                nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
            ispec_save = 0
            ## the output of find_specular_points.c does not start at the same time as the propagation 
            j = 0
            while ( datetime.strptime(read_file_specular[nb_lines_header_output_file_spec].split()[0], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(self.time_sat[j], "%Y-%m-%dT%H:%M:%S") ):
                j = j +1
                self.time_spec_sublist.append('')
                self.name_spec_between_list_and_sublist.append('')
            j = j-1
            while (ispec_save < len(read_file_specular)-1-self.nb_spec_pts):
                j = j + 1
                self.time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
                self.time_spec_sublist_temp_ini = datetime.strptime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
                self.time_spec_sublist.append(datetime.strftime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
                self.name_spec_sublist = []
                self.lon_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1]) #CHANGED
                if self.lon_spec[0,i,j] > 180:
                    self.lon_spec[0,i,j] = self.lon_spec[0,i,j] - 360.
                self.lat_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2]) #CHANGED
                self.gain_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3]) #CHANGED
                self.name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[4]) #CHANGED
                ispec = 1
                while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == self.time_spec_sublist_temp_ini):
                    self.lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
                    if self.lon_spec[ispec,i,j] > 180:
                        self.lon_spec[ispec,i,j] = self.lon_spec[ispec,i,j] - 360.
                    self.lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
                    self.gain_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
                    self.name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4])
                    ispec = ispec + 1
                    if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular) - 1):
                        break
                ispec_save = ispec + ispec_save
                self.name_spec_between_list_and_sublist.append(self.name_spec_sublist)
            # if up to here we still ahve read the entire spec file
            if ( datetime.strptime(self.time_spec_sublist[-1], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(read_file_specular[len(read_file_specular)-2].split()[0], "%Y-%m-%dT%H:%M:%S") ):
                j = j + 1
                first_spec_of_last_time_step = +ispec_save
                self.time_spec_sublist_temp_ini = read_file_specular[first_spec_of_last_time_step].split()[0]
                self.time_spec_sublist_temp_ini = datetime.strptime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
                self.time_spec_sublist.append(datetime.strftime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
                self.name_spec_sublist = []
                self.lon_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[1])
                if self.lon_spec[0,i,j] > 180:
                    self.lon_spec[0,i,j] = self.lon_spec[0,i,j] - 360.
                self.lat_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[2])
                self.gain_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[3])
                ispec = 1
                while (datetime.strptime(read_file_specular[first_spec_of_last_time_step+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == self.time_spec_sublist_temp_ini):
                    self.lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
                    if self.lon_spec[ispec,i,j] > 180:
                        self.lon_spec[ispec,i,j] = self.lon_spec[ispec,i,j] - 360.
                    self.lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
                    self.gain_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
                    if (first_spec_of_last_time_step+ispec < len(read_file_specular) - 1):
                        ispec = ispec + 1
                    else: 
                        break
                self.name_spec_between_list_and_sublist.append(self.name_spec_sublist)

            ## the output of find_specular_points.c does not end at the same time as the propagation 
            j_end = 0
            while ( datetime.strptime(self.time_spec_sublist[-1-j_end], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(self.time_sat[-1-j_end], "%Y-%m-%dT%H:%M:%S") ):
                j_end = j_end +1
                self.time_spec_sublist.append('')
                self.name_spec_between_list_and_sublist.append('')
            self.time_spec.append(self.time_spec_sublist)
        # Build the tuples for the visualization of the specular points
            for k in range(self.nb_spec_pts):
                self.specular = namedtuple('specular',('name',) +  self.point._fields  + self.color._fields + ('point_plot',))
                self.specular_list.append(self.specular)
                # initial position                    
                self.specular_list[k+i*self.nb_spec_pts].x, self.specular_list[k+i*self.nb_spec_pts].y =  m(self.lon_spec[k,i,self.init_index], self.lat_spec[k,i,self.init_index])
                self.specular_list[k+i*self.nb_spec_pts].point_plot = m.plot([],[], marker='o', markersize = 2, color = self.satColors[i], fillstyle = 'none', mew = 2)[0]

            self.name_spec.append(self.name_spec_between_list_and_sublist)


        ###Calculate Default Ground traces###
        self.calcTraces(0, self.nb_steps_interpolation)

        ###Build Ground Stations
        self.gs_name = ['Hawaii', 'Chile', 'Austrailia'] #Note: Order is important for indexing purposes
        self.gs_lat = [21, -33, -31]
        self.gs_lon = [-157, -70, 115]
        self.gsRefs = {}
        self.nb_gs = len(self.gs_name)

        ###Build Storms
        self.forecastTimes = []
        self.forecastLats = []
        self.forecastLons = []
        self.forecastColors = []
        self.forecastRadii = []
        
        self.forecastTimeIdxs = []
        self.forecastTimestamps = []
        self.mapTextRefs = [] #For use in plot storms, probably shouldn't be declared here
        self.activeStorms = []

        self.build_raw_storms()
# v1.2
#        self.generate_interpolators()
        self.interpolate_storms()

        self.forecastArtistlist = [[] for x in range(self.nb_storms)]
        self.trackerArtistlist = [[] for x in range(self.nb_storms)]
        self.trajArtistlist = []

        ###Build Map Time display###
        self.mapTextCorrectionFactor = 2
        #CurTime clock + Display start stop location parameters
        self.pos_display_text_ax_map_x = self.max_lon[self.inZoom] - self.mapTextCorrectionFactor
        self.pos_display_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.display_text_ax_map = self.ax_map.text(self.pos_display_text_ax_map_x, self.pos_display_text_ax_map_y,'',horizontalalignment ='right', weight = 'bold',fontsize = 12, color = 'r')
        #Warning time (forecast time) display
        self.pos_warning_text_ax_map_x = self.min_lon[self.inZoom] + self.mapTextCorrectionFactor
        self.pos_warning_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.time_warning_text_ax_map = self.ax_map.text(self.pos_warning_text_ax_map_x, self.pos_warning_text_ax_map_y,'',horizontalalignment ='left', weight = 'bold',fontsize = 12, color = 'r')
        

        ###Build Sat/GS interaction data structures###
        self.calcGSinteraction()

        ###Build Spec/Storm interaction data structures###
        #NOTE:Used for generating the spec/storm interaction report, but is currently unavailable due to time required
        #to generate report(kcarlsen - 4/30/17)
        #self.findSpecStormPasses()

        ###Call initializing function###
        self.init_data(m)

    def build_raw_storms(self):
        #This list holds tuples of raw forecast data
        self.rawForecasts = []
        #This list holds tuples of (stormName, warningTime) for map
        self.stormInfo = []

        #Need start/end times when deciding indices
        startTime = datetime.strptime(self.time_sat[0], "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.time_sat[-1], "%Y-%m-%dT%H:%M:%S")

        stormFiles = []
        #Get all forecast filenames in ../input/storm_forecasts directory
        if os.path.isdir(self.storm_dir):
            stormFiles = [f for f in os.listdir(self.storm_dir) if os.path.isfile(os.path.join(self.storm_dir, f))]
        else:
            logging.warning("Could not find storm directoy %s\n" \
                + "This may occur if there was an error when downloading storms or if the directory structure" \
                + " for this analysis interval has been altered.", self.storm_dir)

        for file in range(len(stormFiles)):
            STORM_FILEPATH = self.storm_dir + '/' + stormFiles[file]
            stormFile = open(STORM_FILEPATH, 'r')
            #print self.storm_dir + '/' + stormFiles[file]
            forecastLines = stormFile.readlines()

            if len(forecastLines) > 1: #Don't include single forecasts
                #tmpForecasts allows us to append on a per-storm basis
                tmpForecasts = []
                #tmpTimes allows us to only look at 34 kt radii for each storm
                tmpTimes = []
                #tmpLons lets us address duplicate longitudes
                tmpLons = []
                if(STORM_FILEPATH[-3:] == 'txt' and os.stat(STORM_FILEPATH).st_size != 0):
                    #Need first and last cast to determine if storm is in interval
                    firstCast = [x.strip() for x in forecastLines[0].split(',')]
                    firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                    lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                    lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[3]))

                    #TEMP FIX - Checks if any forecast for storm falls in analysis interval
                    #TODO: Fully iron out best way to weed out storms and reduce code duplication here
                    inIntervalFlag = 0
                    for i in range(len(forecastLines)):
                        cast = [x.strip() for x in forecastLines[i].split(',')]
                        castTime = datetime.strptime(cast[2], "%Y%m%d%H") + timedelta(hours=int(cast[3]))
                        if castTime >= startTime and castTime < endTime:
                            inIntervalFlag = 1
                            break #Stop once we know at least one forecast is in interval

                    if inIntervalFlag:
                        logging.debug("Storm %s in analysis interval", firstCast[0])
                        #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                        self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                        #For each line in forecast file
                        for iStorm in range(len(forecastLines)):
                            noaa = [x.strip() for x in forecastLines[iStorm].split(',')]
                            #Build basic storm parameters
                            if noaa[3] not in tmpTimes: #So we only read forecast one time 
                                #Timestamp
                                ts = datetime.strptime(noaa[2], "%Y%m%d%H") + timedelta(hours=int(noaa[3])) 

                                #Index
                                if ts < startTime:
                                    sIdx = -2 #Arbitrary flag 
                                elif ts > endTime:
                                    sIdx = -1 #Arbitrary flag
                                else:
                                    sIdx = int(self.time_sat.index(ts.isoformat()))

                                #Lat
                                lat = float(noaa[4])
                                
                                #Lon
                                lon = float(noaa[5])

                                #Must be no duplicates to satisfy spline condition
                                while lon in tmpLons:
                                    logging.debug("Found %s in tmp lons, altering.", str(lon))
                                    lon -= 0.3
                                tmpLons.append(lon)

                                #Storm type
                                typ = noaa[9]

                                #Size
                                windRadii = [int(noaa[12]), int(noaa[13]), int(noaa[14]), int(noaa[15])]
                                size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding
                                
                                tmpTimes.append(noaa[3]) 

                                newForecast = (ts, sIdx, lat, lon, typ, size)
                                
                                tmpForecasts.append(newForecast)
                        if tmpForecasts:
                            #Create list of lists of tuples (one list for each active storm)
                            self.rawForecasts.append(tmpForecasts)
                
                elif(STORM_FILEPATH[-3:] == 'fst' and os.stat(STORM_FILEPATH).st_size != 0):
                    firstCast = [x.strip() for x in forecastLines[0].split(',')]
                    firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                    lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                    lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[5]))
                    
                    inIntervalFlag = 0
                    for i in range(len(forecastLines)):
                        cast = [x.strip() for x in forecastLines[i].split(',')]
                        castTime = datetime.strptime(cast[2], "%Y%m%d%H") + timedelta(hours=int(cast[3]))
                        if castTime >= startTime and castTime < endTime:
                            inIntervalFlag = 1
                            break #Stop once we know at least one forecast is in interval

                    #Check if storm falls at all in display interval
                    #TEMP FIX
                    #TODO: Fully iron out best way to weed out storms and reduce code duplication here
                    if inIntervalFlag:
                        #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                        self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                        #For each line in forecast file
                        for iStorm in range(len(forecastLines)):
                            jtwc = [x.strip() for x in forecastLines[iStorm].split(',')]
                            #Build basic storm parameters
                            if jtwc[5] not in tmpTimes: #So we only read forecast one time 
                                #Timestamp
                                ts = datetime.strptime(jtwc[2], "%Y%m%d%H") + timedelta(hours=int(jtwc[5])) 

                                #Index
                                if ts < startTime:
                                    sIdx = -2 #Arbitrary flag 
                                elif ts > endTime:
                                    sIdx = -1 #Arbitrary flag
                                else:
                                    sIdx = int(self.time_sat.index(ts.isoformat()))

                                #Lat
                                lat = jtwc[6]
                                if lat[-1:] == 'N':
                                    lat = lat.replace("N","")
                                    lat = float(lat)/10.0
                                else:
                                    lat = lat.replace("S","")
                                    lat = -1*float(lat)/10.0
                                

                                #Lon
                                lon = jtwc[7]
                                if lon[-1:] =='E':
                                    lon = lon.replace("E","")
                                    lon = float(lon)/10.0
                                else:
                                    lon = lon.replace("W","")
                                    lon = -1*float(lon)/10.0

                               #Must be no duplicates to satisfy spline condition
                                while lon in tmpLons:
                                    logging.debug("Found %s in tmp lons, altering.", str(lon))
                                    lon -= 0.3
                                tmpLons.append(lon)

                                #Storm type
                                typ = jtwc[10]

                                #Size
                                windRadii = [int(jtwc[13]), int(jtwc[14]), int(jtwc[15]), int(jtwc[16])]
                                size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding
                                
                                tmpTimes.append(jtwc[5]) 

                                newForecast = (ts, sIdx, lat, lon, typ, size)
                                
                                tmpForecasts.append(newForecast)
                        if tmpForecasts:
                            #Create list of lists of tuples (one list for each active storm)
                            self.rawForecasts.append(tmpForecasts)

        self.nb_storms = len(self.rawForecasts)
#V1.2

    def do_interpolate( self, system ):
        from scipy import interpolate
        lats = [x[2] for x in system]
        lons = [x[3] for x in system]
        if len( lats ) == 3:
            tck,u=interpolate.splprep([lons,lats],s=0.0, k=2)
        elif len( lats ) == 2:
            tck,u=interpolate.splprep([lons,lats],s=0.0, k=1)
        else:
            tck,u=interpolate.splprep([lons,lats],s=0.0)
        # 40000 > num values rounded to tenths place between -180 and 180
        newLons,newLats = interpolate.splev(np.linspace(0,1,40000),tck) 
        d = {}
        for i,newLon in enumerate(newLons):
            d[newLon] = newLats[i]
        return( d )
        
    def interpolate_storms(self):
        finalTimestamps = [[] for x in range(self.nb_storms)]
        finalIdxs = [[] for x in range(self.nb_storms)]
        finalLats = [[] for x in range(self.nb_storms)]
        finalLons = [[] for x in range(self.nb_storms)]
        finalTypes = [[] for x in range(self.nb_storms)]
        finalSizes = [[] for x in range(self.nb_storms)]
        
        timeSkip = 720 #12 minutes in seconds
        startTime = datetime.strptime(self.time_sat[0], "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.time_sat[-1], "%Y-%m-%dT%H:%M:%S")

        systemIdx = 0
        for system in self.rawForecasts:
            if len(system) > 1:
                for i in range(len(system) - 1):
                    #Figure out number of steps in this interval
                    intervals = int(divmod((system[i+1][0] - system[i][0]).total_seconds(), timeSkip)[0])
                    for j in range(intervals):
                        #Build new timestamps
                        fts = system[i][0] + j*timedelta(seconds=timeSkip)
                        finalTimestamps[systemIdx].append(fts)

                        #Retrieve new indices
                        if fts < startTime:
                            finalIdxs[systemIdx].append(-2) #Arbitrary flag 
                        elif fts > endTime:
                            finalIdxs[systemIdx].append(-1) #Arbitrary flag
                        else:
                            finalIdxs[systemIdx].append(int(self.time_sat.index(fts.isoformat())))

                    #Generate new lons
                    newLons = np.linspace(system[i][3], system[i+1][3], intervals, endpoint=False)
                    finalLons[systemIdx].extend(newLons)

                    #Generate new types
                    finalTypes[systemIdx].extend([system[i][4]] * intervals)

                    #Generate new sizes
                    finalSizes[systemIdx].extend(np.linspace(system[i][5], system[i+1][5], intervals, endpoint=False))
            
            #Generate new lats from interpolator
            finalLats[systemIdx] = []
            dLatLons = self.do_interpolate( system )
            for finalLon in finalLons[systemIdx]:
                lookup = dLatLons.get( finalLon, -9999 )
                if lookup == -9999:
                    arr = np.array(dLatLons.keys())
                    closest = arr[ np.abs( arr-finalLon ).argmin() ]
                    finalLats[systemIdx].append( dLatLons[closest] )
                else:
                    finalLats[systemIdx].append( dLatLons[finalLon] )
 
# Adding code to sort finalLats/Lons by euclidian distance so storms update smoothly
            combined = zip( finalLons[systemIdx],finalLats[systemIdx] )
            whatsLeft = combined[:]
            myPt = whatsLeft[0]
            final = []
            final.append( myPt)
            whatsLeft = [ i for i in whatsLeft if i!=myPt ]
            while whatsLeft:
                nextPt = whatsLeft[ distance.cdist( [myPt],np.array(whatsLeft) ).argmin()]
                final.append( nextPt )
                whatsLeft = [ i for i in whatsLeft if i!=nextPt ]
                myPt = nextPt
            newLons,newLats = zip(*final)
            finalLons[systemIdx] = list(newLons)
            finalLats[systemIdx] = list(newLats)
                
                        
            systemIdx += 1
       
       
        #Find start/end indices for each storm
        self.stormStart = []
        self.stormEnd = []
        self.curStorm = []

        for i in range(self.nb_storms):
            self.stormStart.append(bisect(finalIdxs[i], -2)) #Find first value not before sim start
            self.curStorm.append(self.stormStart[i])
            try:
                self.stormEnd.append(finalIdxs[i].index(-1)) #Find first value past end of interval
            except ValueError:
                self.stormEnd.append(len(finalIdxs[i]) - 1)

        self.finalForecasts = []

        for i in range(self.nb_storms):
            self.finalForecasts.append(zip(finalTimestamps[i], finalIdxs[i], finalLats[i], finalLons[i], finalTypes[i], finalSizes[i]))


    ### Data/Simulation managing functions ###
    def init_data(self, m):
        #SATELLITES AND SPECULAR POINTS
        self.tuple_spacecraft_point_to_plot = ()
        self.tuple_specular_point_to_plot = ()
        for isat in range(self.nb_satellites):
            self.spacecraft_list[isat].point_plot.set_data([], []) # set_data must be applied to a line object
            self.tuple_spacecraft_point_to_plot = self.tuple_spacecraft_point_to_plot + (self.spacecraft_list[isat].point_plot,)
            for ispec in range(self.nb_spec_pts):
                self.specular_list[ispec+isat*self.nb_spec_pts].point_plot.set_data([], [])
                self.tuple_specular_point_to_plot = self.tuple_specular_point_to_plot + (self.specular_list[ispec+isat*self.nb_spec_pts].point_plot,)

        #Prepare for animation
        self.firstTime = self.time_sat[0] #First and last time maintained separately to maintain absolute bounds
        self.lastTime = self.time_sat[-1]

        #These values are strings
        self.start_visu = self.firstTime
        self.stop_visu = self.lastTime

        #These values are indices
        self.when_start_visu = self.time_sat.index(self.start_visu)
        self.when_stop_visu = self.time_sat.index(self.stop_visu)
        self.curSample = self.when_start_visu

        self.tuple_spacecraft_and_specular = self.tuple_spacecraft_point_to_plot + self.tuple_specular_point_to_plot
        
        logging.info("Analysis Interval: %s  -> %s", self.time_sat[self.when_start_visu], self.time_sat[self.when_stop_visu])

        #TIME
        self.display_text_ax_map.set_text('')
        self.time_warning_text_ax_map.set_text("Reported Forecast Times \n --- \n" + formatWarningTime(self.stormInfo))
        self.tuple_spec_text_ax_map = ()

        #GS_INTERACTION
        self.GSinteractionIdx = [0] * self.nb_satellites #We need 8 indices to mark our progress, one for each sat

    def advSim(self, plotMap):
        self.advSimtuple_spacecraft_point_to_plot = ()
        self.advSimtuple_specular_point_to_plot = ()
        
        #This is where the sats are updated every iteration during animation
        #To not show, only run inner code block for satellites marked as 'true'
        for isat in range(self.nb_satellites):
            # SPACECRAFT
            lon_a, lat_a = self.lon_sat[isat,self.curSample], self.lat_sat[isat,self.curSample]
            self.spacecraft_list[isat].x,  self.spacecraft_list[isat].y = plotMap(lon_a, lat_a)
            self.spacecraft_list[isat].point_plot.set_data(self.spacecraft_list[isat].x, self.spacecraft_list[isat].y)       
            # SPECULAR POINTS
            kspec = 0
            # if there is less specular points than self.nb_spec_pts for for this iteration and satellite, then the old specular point stays on the image. So for loop blow is to make sure that does not happen (there are smarter ways to do that...)
            for ispec in range(self.nb_spec_pts):
                self.specular_list[ispec + isat*self.nb_spec_pts].point_plot.set_data([], []) 
            for ispec in range(len(self.name_spec[isat][self.curSample])):
                lon_spec_a, lat_spec_a = self.lon_spec[ispec,isat,self.curSample], self.lat_spec[ispec,isat,self.curSample]
                self.specular_list[ispec + isat*self.nb_spec_pts].x,  self.specular_list[ispec + isat*self.nb_spec_pts].y = plotMap(lon_spec_a, lat_spec_a)
                self.specular_list[ispec+isat*self.nb_spec_pts].point_plot.set_data(self.specular_list[ispec+isat*self.nb_spec_pts].x, self.specular_list[ispec+isat*self.nb_spec_pts].y)
                self.advSimtuple_specular_point_to_plot = self.advSimtuple_specular_point_to_plot + (self.specular_list[ispec+isat*self.nb_spec_pts].point_plot,)

            # GS INTERACTION
            if self.animDirection == 1: #Animation in forward
                if self.GSinteractionIdx[isat] < len(self.overGSstart[isat]) and datetime.strptime(self.time_sat[self.curSample], '%Y-%m-%dT%H:%M:%S') > self.overGSstart[isat][self.GSinteractionIdx[isat]]:
                    self.GSinteractionIdx[isat] += 1
            elif self.animDirection == -1: #Animation in reverse
                if (self.GSinteractionIdx[isat] > 0 and datetime.strptime(self.time_sat[self.curSample], '%Y-%m-%dT%H:%M:%S') < self.overGSstart[isat][self.GSinteractionIdx[isat] - 1]):
                    logging.debug("Decrement for %s", str(isat))
                    self.GSinteractionIdx[isat] -= 1


        #TIME
        self.display_text_ax_map.set_text("Start: " + formatDisplayTime(self.time_sat[self.when_start_visu]) + "\n"  +  "Current: "  + formatDisplayTime(self.time_sat[self.curSample]) + "\n" + "End:   " + formatDisplayTime(self.time_sat[self.when_stop_visu]))

        #STORMS
        if self.stormsShowed:
            for storm in range(self.nb_storms):            
                if self.curStorm[storm] != self.stormEnd[storm] and self.curSample >= self.finalForecasts[storm][self.curStorm[storm] + 1][1]:
                    self.curStorm[storm] += 1
                    self.removeFirstStorm(storm)
                    self.updateFirstStorm(storm)
 
    def incrementIndex(self):

        if(self.curSample == self.when_stop_visu):
            self.curSample = self.when_start_visu
            self.forwardLoopReset()

        if(self.curSample + self.speedFactor < self.when_stop_visu):
            self.curSample += self.speedFactor
        else:
            self.curSample = self.when_stop_visu

    def decrementIndex(self):
        if(self.curSample == self.when_start_visu):
            self.curSample = self.when_stop_visu
            self.findGSinteractionIdx()

        if(self.curSample - self.speedFactor > self.when_start_visu):
            self.curSample -= self.speedFactor
        else:
            self.curSample = self.when_start_visu

    def updateCurStorm(self):
        self.validFlags = []
        for storm in range(self.nb_storms):
            idxs = [x[1] for x in self.finalForecasts[storm]]
            logging.debug("idxs: %s", str(idxs))

            stormIdx, validFlag = self.find_le(idxs, self.curSample)
            self.validFlags.append(validFlag)
            self.curStorm[storm] = idxs.index(stormIdx)
            logging.debug("curStorm = %s",self.curStorm[storm])

    #TODO: Create a reverseLoopReset()

    def forwardLoopReset(self):
        self.findGSinteractionIdx()
        if self.stormsShowed:
            self.removeStorms()
            # self.updateCurStorm()
            self.plotStorms(self.plotMap)


    ### Map Managing functions ###
    def radius_for_tissot(self,dist_km):
        return np.rad2deg(dist_km/6378.) # this is calculating using the Haversine formula
        
    def rebuildTimeDisplay(self, zoom):

        self.inZoom = zoom
        #Remove old text
        self.display_text_ax_map.remove()
        self.time_warning_text_ax_map.remove()

        self.pos_display_text_ax_map_x = self.max_lon[self.inZoom] - self.mapTextCorrectionFactor
        self.pos_display_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.display_text_ax_map = self.ax_map.text(self.pos_display_text_ax_map_x, self.pos_display_text_ax_map_y,'',horizontalalignment ='right', weight = 'bold',fontsize = 15, color = 'r')

        #Interval display
        self.pos_warning_text_ax_map_x = self.min_lon[self.inZoom] + self.mapTextCorrectionFactor
        self.pos_warning_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.time_warning_text_ax_map = self.ax_map.text(self.pos_warning_text_ax_map_x, self.pos_warning_text_ax_map_y,'',horizontalalignment ='left', weight = 'bold',fontsize = 15, color = 'r')
        #Set warning text here since it doesn't change over the course of the animation
        self.time_warning_text_ax_map.set_text("Reported Warning Times \n ---\n" + formatWarningTime(self.stormInfo))

    def plotGroundStation(self, plotMap, idx):
        #self.gsRefs[self.gs_name[idx]] = plotMap.plot(self.gs_lon[idx], self.gs_lat[idx], color='grey', marker='o', markersize=100, fillstyle='none', label=self.gs_name[idx],linestyle=':')
        #self.gsRefs[self.gs_name[idx]] = plotMap.tissot(self.gs_lon[idx], self.gs_lat[idx], self.radius_for_tissot(2000), 256, facecolor='none', linestyle='dashdot', linewidth=3)
        self.gsRefs[self.gs_name[idx]] = plotMap.tissot(self.gs_lon[idx], self.gs_lat[idx], self.radius_for_tissot(2000), 100, linestyle='dashdot', fill=False, edgecolor='grey')
         
    def removeGroundStation(self, idx):
        doomedGS = self.gsRefs.pop(self.gs_name[idx]) #This pop refers to the dictionary holding map objects
        doomedGS.remove() #This pop refers to the artist managing the map object

    def find_le(self, inList, x):
        #'Find rightmost value less than or equal to x'
        i = bisect_right(inList, x)
        if i and inList[i-1] >=0:
            return inList[i-1], 1
        else:
            return inList[0], 0

    def plotStorms(self, plotMap):
        self.activeStorms = []
        #This save and then repass in forwardloopreset seems unnecessary - improve?
        self.plotMap = plotMap
        first_in_past = True

        self.updateCurStorm()
        
        for storm in range(self.nb_storms):
            if self.validFlags[storm]:
                for idx in range(self.curStorm[storm], len(self.finalForecasts[storm])):
                    self.forecastArtistlist[storm].append(plotMap.tissot(self.finalForecasts[storm][idx][3], self.finalForecasts[storm][idx][2], radius_for_tissot(self.finalForecasts[storm][idx][5]), 256, facecolor=self.colorByIntensity_faded[self.finalForecasts[storm][idx][4]], edgecolor='none'))
                
                self.updateFirstStorm(storm)
            else:
                for idx in range(len(self.finalForecasts[storm])):
                        self.forecastArtistlist[storm].append(plotMap.tissot(self.finalForecasts[storm][idx][3], self.finalForecasts[storm][idx][2], radius_for_tissot(self.finalForecasts[storm][idx][5]), 256, facecolor=self.colorByIntensity_faded[self.finalForecasts[storm][idx][4]], edgecolor='none'))
        
        self.stormsShowed = True
             

    def removeStorms(self):
        for storm in range(self.nb_storms):
            while len(self.forecastArtistlist[storm]) > 0:
                self.forecastArtistlist[storm][len(self.forecastArtistlist[storm]) - 1].remove()
                self.forecastArtistlist[storm].pop()
            self.stormsShowed = False

    def removeFirstStorm(self, storm):
        if len(self.forecastArtistlist[storm]) > 0:
            self.forecastArtistlist[storm][0].remove()
            self.forecastArtistlist[storm].pop(0)

    def updateFirstStorm(self, storm):
        if len(self.forecastArtistlist[storm]) > 0:
            self.forecastArtistlist[storm][0].set_zorder(2)
            self.forecastArtistlist[storm][0].set_facecolor(self.colorByIntensity[self.finalForecasts[storm][self.curStorm[storm]][4]])
            self.forecastArtistlist[storm][0].set_edgecolor('black')

    def plotTrajectory(self, plotMap):
        for storm in range(self.nb_storms):
            stormLons = [x[3] for x in self.finalForecasts[storm]]
            stormLats = [x[2] for x in self.finalForecasts[storm]]
            self.trajArtistlist.append(plotMap.scatter(stormLons, stormLats, zorder=5, s = 1, facecolor='0.5', lw = 0, alpha=0.5))

        
    #def monotonic(self,x):
    #    dx = np.diff(x)
    #    return np.all(dx <= 0) or np.all(dx >= 0)
    #    
    #def plotTrajectory(self, plotMap):
    #    for storm in range(self.nb_storms):
    #        stormLons = [x[3] for x in self.finalForecasts[storm]]
    #        stormLats = [x[2] for x in self.finalForecasts[storm]]
    #        
    #        combined = zip( stormLons,stormLats )
    #        allL = []
    #        thisLonList = []
    #        thisLonLatList = []
    #        i=0
    #        while i < len(combined):
    #            mytup = combined[i]
    #            myLon = mytup[0]
    #            myLat = mytup[1]
    #            
    #            thisLonList.append(myLon)
    #            thisLonLatList.append( mytup )
    #            
    #            if self.monotonic(thisLonList):
    #                i+=1
    #            else:
    #                allL.append(thisLonLatList[:-1])
    #                thisLonList = []
    #                thisLonLatList = []
    #        allL.append( thisLonLatList )
    #        for l in allL[0:2]:
    #            myLons,myLats = zip(*l)
    #            myLons = list(myLons)
    #            myLats = list(myLats)
    #            self.trajArtistlist.append(plotMap.scatter(myLons, myLats,zorder=5))
        
    def removeTrajectory(self):
        logging.debug("%d", self.nb_storms)
        for storm in range(self.nb_storms):
            doomedTraj = self.trajArtistlist.pop(0)
#            doomedTraj.pop(0).remove()
            doomedTraj.remove()

    #Unused in first revision of SIFT(kcarlsen - 4/30)

    # def plotForecastText(self):
    #     for idx in self.activeStorms:
    #         self.mapTextRefs.append(self.ax_map.text(self.forecastLons[idx] + 5, self.forecastLats[idx] + 8, '', rotation=45, horizontalalignment ='right', fontsize = 8, color = 'black').set_text(self.forecastTimestamps[idx]))

    # def removeForecastText(self):
    #     while len(self.mapTextRefs) > 0:
    #         self.mapTextRefs[len(self.mapTextRefs) - 1].remove(0)
    #         self.mapTextRefs.pop()

    def calcTraces(self, startSample, endSample):
        self.satTraceX = [[] for i in range(self.nb_satellites)]
        self.satTraceY = [[] for i in range(self.nb_satellites)]
        self.specTraceX = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceY = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceGain = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]

        for i in range(self.nb_satellites):
            for j in range(startSample, endSample):
                self.satTraceX[i].append(self.lon_sat[i, j])
                self.satTraceY[i].append(self.lat_sat[i, j])

            for k in range(self.nb_spec_pts):
                for step in range(startSample, endSample-1): #Not sure about this -1, but otherwise there is an extra 0
                    self.specTraceX[i][k].append(self.lon_spec[k,i,step])
                    self.specTraceY[i][k].append(self.lat_spec[k,i,step])
                    if self.gain_spec[k,i,step] < 5:
                        self.specTraceGain[i][k].append(0.5)
                    elif self.gain_spec[k,i,step] < 10:
                        self.specTraceGain[i][k].append(5)
                    else:
                        self.specTraceGain[i][k].append(15)

    def plotSatTrace(self, plotMap, satNum):
        self.satTraceRefs[self.satNames[satNum]] = plotMap.plot(self.satTraceX[satNum], self.satTraceY[satNum], 'red', linestyle='', marker='.', markersize=1)

    def removeSatTrace(self, satNum):
        doomedTrace = self.satTraceRefs.pop(self.satNames[satNum])
        doomedTrace.pop(0).remove()

    def plotSpecTrace(self, plotMap, satNum):
        self.specTraceRefs[self.satNames[satNum]] = []
        for idx in range(len(self.specNames)):
            self.specTraceRefs[self.satNames[satNum]].append(plotMap.plot(self.specTraceX[satNum][idx], self.specTraceY[satNum][idx], color=self.satColors[satNum], linestyle='', marker='.', markersize=1))

    def removeSpecTrace(self, satNum):
        doomedTraces = self.specTraceRefs[self.satNames[satNum]]
        for trace in doomedTraces:
            trace.pop(0).remove()

    def plotSpecGainTrace(self, plotMap, satNum):
        self.specGainTraceRefs[self.satNames[satNum]] = []
        for idx in range(len(self.specNames)):
            #s is the scale factor for scatter plot
            #c is color
            #zorder=3 allows gain lines to plot on top of storms
            self.specGainTraceRefs[self.satNames[satNum]].append(plotMap.scatter(self.specTraceX[satNum][idx], self.specTraceY[satNum][idx], s=self.specTraceGain[satNum][idx], color=self.satColors[satNum], zorder=3)) #This allows variable markersize but is slow

    def removeSpecGainTrace(self, satNum):
        doomedTraces = self.specGainTraceRefs[self.satNames[satNum]]
        for trace in doomedTraces:
            trace.remove()


    def calcGSinteraction(self):
        self.overGSstart = []
        self.overGSend = []
        self.overGSdur = []
        self.overGSwhich = []

        for i in range(self.nb_satellites):
            #Read in file
            data_filename = self.gsInteraction_filepath + self.output_filename_propagator[i]
            data_file = open(data_filename, 'r')
            data_lines = data_file.readlines()
            #Pop off Header and blank line
            del data_lines[0]
            del data_lines[0]

            temp_start = []
            temp_end = []
            temp_dur = []
            temp_which = []

            for j in range(len(data_lines)):
                line = data_lines[j].split()

                start_time = line[1] + "T" + line[2]
                temp_start.append(datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S'))

                end_time = line[4] + "T" + line[5]
                temp_end.append(datetime.strptime(end_time, '%Y-%m-%dT%H:%M:%S'))

                temp_dur.append(temp_end[j] - temp_start[j]) #Appends as timedelta object

                temp_which.append(line[12])

            self.overGSstart.append(temp_start)
            self.overGSend.append(temp_end)
            self.overGSdur.append(temp_dur)
            self.overGSwhich.append(temp_which)

    def findSpecStormPasses(self):
        self.storm_names = []
        self.specPassStart = []
        self.specPassGps = []
        self.specPassDist = []
        self.specPassIn = []
        
        for i in range(self.nb_satellites):
            data_filename = self.specInteraction_filepath + self.output_filename_propagator[i]
            data_file = open(data_filename, 'r')
            data_lines = data_file.readlines()
            #Pop off Header line
            del data_lines[0]
            
            if len(self.storm_names) == 0:
                line = data_lines[0].split()
                for m in range(self.nb_storms):
                    name = line[5+2*m].split('_')[1]
                    self.storm_names.append(name)

            del data_lines[0]
            del data_lines[0]

            temp_start = [[] for x in range(self.nb_storms)]
            temp_gps = [[] for x in range(self.nb_storms)]
            temp_dist = [[] for x in range(self.nb_storms)]
            temp_in = [[] for x in range(self.nb_storms)]

            for j in range(len(data_lines) - 1):
                line = data_lines[j].split()
                
                if datetime.strptime(line[0], '%Y-%m-%dT%H:%M:%S') > datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S') and datetime.strptime(line[0], '%Y-%m-%dT%H:%M:%S') < datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S'):            
                    for k in range(self.nb_storms):
                        if line[6+2*k] == "1": #Checks that spec is in storm
                            temp_start[k].append(line[0])
                            temp_gps[k].append(line[4])
                            temp_dist[k].append(line[5+2*k])
                            temp_in[k].append(line[6+2*k])

            self.specPassStart.append(temp_start)
            self.specPassGps.append(temp_gps)
            self.specPassDist.append(temp_dist)
            self.specPassIn.append(temp_in)

# v1.2

    def generateOverpassReport(self, lowerLon, upperLon, lowerLat, upperLat, start, end, includePoints):

        def xformLon( lon ):
            if lon <= 180:
                return( lon )
            elif lon > 180:
                return (lon-360)
        
        def xformGain( gain ):
            gain = float(gain)
            if gain < 5.0:
                return('LOW')
            elif gain >= 5.0 and gain <=10:
                return('MEDIUM')
            elif gain > 10:
                return('HIGH')
            else:
                return(gain)
                
        dNav = {'NAVSTAR_10': 'N/A',
        'NAVSTAR_11': 'N/A',
        'NAVSTAR_13': 'N/A',
        'NAVSTAR_14': 'N/A',
        'NAVSTAR_15': 'N/A',
        'NAVSTAR_16': 'N/A',
        'NAVSTAR_17': 'N/A',
        'NAVSTAR_18': 'N/A',
        'NAVSTAR_19': 'N/A',
        'NAVSTAR_20': 'N/A',
        'NAVSTAR_21': 'N/A',
        'NAVSTAR_22': 'N/A',
        'NAVSTAR_23': 'N/A',
        'NAVSTAR_24': 'N/A',
        'NAVSTAR_25': 'N/A',
        'NAVSTAR_26': 'N/A',
        'NAVSTAR_27': 'N/A',
        'NAVSTAR_28': 'N/A',
        'NAVSTAR_29': 'N/A',
        'NAVSTAR_30': 'N/A',
        'NAVSTAR_31': 'N/A',
        'NAVSTAR_32': 'N/A',
        'NAVSTAR_33': 'N/A',
        'NAVSTAR_34': 'N/A',
        'NAVSTAR_35': 'N/A',
        'NAVSTAR_36': 'N/A',
        'NAVSTAR_37': 'N/A',
        'NAVSTAR_38': 'N/A',
        'NAVSTAR_39': 'N/A',
        'NAVSTAR_43': '13',
        'NAVSTAR_44': 'N/A',
        'NAVSTAR_46': '11',
        'NAVSTAR_47': '20',
        'NAVSTAR_48': '28',
        'NAVSTAR_49': '14',
        'NAVSTAR_50': '18',
        'NAVSTAR_51': '16',
        'NAVSTAR_52': '21',
        'NAVSTAR_53': '22',
        'NAVSTAR_54': '19',
        'NAVSTAR_55': '23',
        'NAVSTAR_56': '2',
        'NAVSTAR_57': '17',
        'NAVSTAR_58': '31',
        'NAVSTAR_59': '12',
        'NAVSTAR_60': '15',
        'NAVSTAR_61': '29',
        'NAVSTAR_62': '7',
        'NAVSTAR_63': 'N/A',
        'NAVSTAR_64': '5',
        'NAVSTAR_65': '25',
        'NAVSTAR_66': '1',
        'NAVSTAR_67': '24',
        'NAVSTAR_68': '27',
        'NAVSTAR_69': '30',
        'NAVSTAR_70': '6',
        'NAVSTAR_71': '9',
        'NAVSTAR_72': '3',
        'NAVSTAR_73': '26',
        'NAVSTAR_74': '8',
        'NAVSTAR_75': '10',
        'NAVSTAR_76': '32',
        'NAVSTAR_9': 'N/A'}

        dSatNum = {}
        dSatNum[1]=5
        dSatNum[2]=4
        dSatNum[3]=2
        dSatNum[4]=1
        dSatNum[5]=8
        dSatNum[6]=6
        dSatNum[7]=7
        dSatNum[8]=3
        
        LLlat = float(lowerLat)
        LLlon = float(lowerLon)
        URlat = float(upperLat)
        URlon = float(upperLon)
        
        startString = start.strftime('%Y-%m-%dT%H-%M-%S')
        endString = end.strftime('%Y-%m-%dT%H-%M-%S')
        startdt = start
        enddt = end

        if not os.path.isdir(self.overpassReportLoc):
            os.makedirs(self.overpassReportLoc)

        file_loc = self.overpassReportLoc + "/" + str(datetime.now().strftime('%Y-%m-%dT%H-%M-%S')) + ".csv"

        if os.path.isfile(file_loc): #Remove old version of report if it exists
            os.remove(file_loc)

        outfile = open(file_loc, 'w')

        files = glob.glob( self.specular_filepath+'*' )

        dfList = []
        for fn in files:
            satNum = os.path.basename(fn).split('.')[0][-1]
            satName = 'CYGFM0'+str(dSatNum[int(satNum)])
            df = pd.read_csv(fn,delim_whitespace=True,skiprows=4,parse_dates=[0])
            df.columns = ['TIME', 'LON_SPEC', 'LAT_SPEC', 'GAIN', 'NAME_GPS']
            df['LON_SPEC'] = df['LON_SPEC'].map( xformLon )
            df['GAIN'] = df['GAIN'].map( xformGain )
            mydf = df.ix[ (df['LAT_SPEC']<= URlat) & (df['LAT_SPEC'] >= LLlat) & (df['LON_SPEC'] <= URlon ) & (df['LON_SPEC'] >= LLlon) ]
            mydf = mydf.ix[ (mydf['TIME'] >= startdt) & (mydf['TIME']<= enddt)  ]
            mydf['CYGFM0#']=str(dSatNum[int(satNum)])
            headFn = glob.glob( self.sat_filepath+'spock_spec_start_*'+str(satNum)+'.txt' )[0]        
            dfHead = pd.read_csv( headFn,delim_whitespace=True,parse_dates=[0])
            dfHead = dfHead[[dfHead.columns[0],dfHead.columns[-1]]]
            dfHead.columns = ['TIME','HEADING']
            dfNew = pd.merge( mydf,dfHead,on='TIME' )
            dfList.append(dfNew)
        df = pd.concat(dfList,axis=0)
        df = df[ ['CYGFM0#','TIME','LAT_SPEC','LON_SPEC','HEADING','GAIN','NAME_GPS'] ]
        df = df.sort_values( by=['CYGFM0#','TIME'] )
        df['NUM_PRN'] = df['NAME_GPS'].map(dNav)
        df['SCIENCE_MODE_Y/N']= 'Y'
        df.to_csv(outfile,index=False)
        
#        os.startfile(os.path.abspath(file_loc))
        logging.info("Overpass report generated.")


    def generateSpecReport(self):
        file_loc = self.specReportLoc + "/" + str(datetime.now().strftime('%Y-%m-%dT%H-%M-%S')) + ".txt"
        outfile = open(file_loc, 'w')
        startTime = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S')
        endTime = datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S')
        self.findSpecStormPasses()

        #Write header
        outfile.write("CYGNSS Specular Flyover Report\n")
        outfile.write("-----\n")
        outfile.write("Generated on: " + str(datetime.now()) + "\n")
        outfile.write("For time range: " + self.start_visu + " to " + self.stop_visu + "\n")
        outfile.write("Number of Satellites in this report: " + str(self.nb_satellites) + "\n")
        outfile.write("-----\n\n")

        for sat in range(self.nb_satellites):
            #outfile.write("###CYGNSS " + str(sat + 1) + "###\n") #Sat+1 converts from 0 to 1 indexing
            outfile.write("###CYG" + cygnss_name_arr[sat] + "###\n")
            for storm in range(len(self.storm_names)):
                outfile.write(self.storm_names[storm])
                outfile.write("\n")
                for spec in range(len(self.specPassStart[sat][storm])):
                    outfile.write(self.specPassStart[sat][storm][spec] + "   " + self.specPassDist[sat][storm][spec] + "   " + self.specPassGps[sat][storm][spec])
                    outfile.write("\n")
                outfile.write("\n")
            outfile.write("---")
            outfile.write("\n")

        outfile.write("-----\n")
        outfile.write("Note: This file was autogenerated by SIFT.\n")

        outfile.close()
        logging.info("Specular/Storm report generated.")

        if sys.platform.startswith('win'):
		os.startfile(file_loc)
        else:
        	if sys.platform == "darwin":
        		opener = "open"        		
        	else:
        		opener = "xdg-open"
        	subprocess.call([opener,file_loc])

    #Sets appropriate gsinteraction to current indices when loop wraps around
    def findGSinteractionIdx(self):
        for i in range(self.nb_satellites):
            try:
                if self.animDirection >= 0:
                    self.GSinteractionIdx[i] = 0 #Reset value to base
                    tmpDT = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S') #start_visu because forward loop
                    while tmpDT > self.overGSstart[i][self.GSinteractionIdx[i]]:
                        self.GSinteractionIdx[i] += 1
                else:
                    self.GSinteractionIdx[i] = len(self.overGSstart[i]) #Reset value to valid max for sat
                    tmpDT = datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S') #end_visu because backward loop
                    while tmpDT < self.overGSstart[i][self.GSinteractionIdx[i]]:
                        logging.debug("In decrement")
                        self.GSinteractionIdx[i] -= 1
                #print "GSinteractionIdx[" + str(i) +"] = " + str(self.GSinteractionIdx[i])
            except IndexError:
                #print "GSinteractionIdx[" + str(i) +"] = " + str(self.GSinteractionIdx[i])
                pass

    def generateGSreport(self):
        if not os.path.isdir(self.gsReportLoc):
            os.makedirs(self.gsReportLoc)

        file_loc = self.gsReportLoc + "/" + str(datetime.now().strftime('%Y-%m-%dT%H-%M-%S')) + ".txt"

        if os.path.isfile(file_loc): #Remove old version of report if it exists
            os.remove(file_loc)

        outfile = open(file_loc, 'w')
        startTime = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S')
        endTime = datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S')

        #Write header
        outfile.write("CYGNSS Ground Station Flyover Report\n")
        outfile.write("-----\n")
        outfile.write("Generated on: " + str(datetime.now()) + "\n")
        outfile.write("For time range: " + self.start_visu + " to " + self.stop_visu + "\n")
        outfile.write("Number of Satellites in this report: " + str(self.nb_satellites) + "\n")
        outfile.write("-----\n")

        #Find first GSInteractionIdx of interest for all sats
        firstGSIdx = self.GSinteractionIdx

        #Find last GSInteractionIdx of intererst for all sats
        lastGSIdx = []
        for i in range(self.nb_satellites):
            tmpIdx = 0
            while tmpIdx < len(self.overGSend[i]) and endTime > self.overGSend[i][tmpIdx]:
                tmpIdx += 1

            lastGSIdx.append(tmpIdx)

        #Write information to file
        for i in range(self.nb_satellites):
            outfile.write(self.satNames[i] + "\n")
            for interactionIdx in range(firstGSIdx[i], lastGSIdx[i]):
                outfile.write("\tGround Station: " + self.overGSwhich[i][interactionIdx] + "  Connection Start: " + str(self.overGSstart[i][interactionIdx]) + 
                    "  Connection End: " + str(self.overGSend[i][interactionIdx]) + "  Duration: " + str(self.overGSdur[i][interactionIdx]) + "\n")
                interactionIdx += 1

            outfile.write("\n")

        outfile.write("-----\n")
        outfile.write("Note: This file was autogenerated by SIFT.  It only includes Ground Station interactions that COMPLETELY fall within the analysis interval given at the top of this file\n")

        outfile.close()

        logging.info("Satellite/Ground Station report generated.")
        if sys.platform.startswith('win'):
        	os.startfile(os.path.abspath(file_loc))
        else:
        	if sys.platform == "darwin":
        		opener = 'open'
        	else:
        		opener =  "xdg-open"
        	subprocess.call([opener,os.path.abspath(file_loc)])




##############################################


class CygnssFrame(wx.Frame):
    """ The main application frame
    """
    title = 'Storm Intersection Forecast Tool'
    zoom = 0
    mapRef = None
    localGSInteractionIdxs = None
    meridians = {}
    parallels = {}
    userZoom_filepath = '../input_sift/user_zoom.txt'
#v1.2
    bullseyes = {}


    def __init__(self):
#v1.2
#        wx.Frame.__init__(self, None, -1, self.title)
#v1.2
        displaySize= wx.DisplaySize()
#        wx.Frame.__init__(self, None, -1, self.title, size = (displaySize[0],displaySize[1]) )
        wx.Frame.__init__(self, None, -1, self.title )
        
#v1.2
        self.initialize_plot_map()
        
        self.satManager = SatelliteManager(self.plotMap, self.max_lon, self.min_lon, self.min_lat, self.height_fig_with_map, self.width_fig_with_map, self.ax_map, self.zoom)
        self.localGSInteractionIdxs = [-1] * self.satManager.nb_satellites
        self.create_main_panel()
        #Icon must be added after creation of panel
        ico = cygIcon.cyg_logo_small.getIcon()
        self.SetIcon(ico)

        self.paused = True
        self.revPaused = True
        self.create_menu()
        

        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)
        self.redraw_timer.Start(50) #Controls how often the UI is redrawn

        #Bind 'X' button to on_exit to ensure proper shut down
        self.Bind(wx.EVT_CLOSE, self.on_exit)
#v1.2
# Add ability to click dots

        cid = self.fig_with_map.canvas.mpl_connect('button_press_event',self.on_press)



    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_exit = menu_file.Append(-1, "&Exit\t Ctrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

        menu_view = wx.Menu()
        zoomMenuItems = []

        #Make menu items
        for zoom in self.zoom_names:
            name = "&" + zoom
            zoomMenuItems.append(menu_view.Append(-1, name, zoom))

        #Bind all items to handler
        for item in zoomMenuItems:
            self.Bind(wx.EVT_MENU, self.on_zoom, item)

        menu_generate = wx.Menu()
        overpassReport = menu_generate.Append(-1, "&Overpass Report", "Overpass Report")
        self.Bind(wx.EVT_MENU, self.on_overpass_report, overpassReport)
        stormReport = menu_generate.Append(-1, "&Specular/Storm Report", "Specular/Storm Report")
        self.Bind(wx.EVT_MENU, self.on_generate_spec_report, stormReport)
        stormReport.Enable(False) #This feature not available in SIFT 1.0, but code is partially written
        gsReport = menu_generate.Append(-1, "&Satellite/Ground Station Report", "Satellite/Ground Station Report")
        self.Bind(wx.EVT_MENU, self.on_generate_gs_report, gsReport)


        menu_help = wx.Menu()
        stormColorKey = menu_help.Append(-1, "&Storm Forecast Color Key", "Storm Forecast Color Key")
        self.Bind(wx.EVT_MENU, self.on_storm_color_key, stormColorKey)
        gainThickKey = menu_help.Append(-1, "&Specular Trace Gain Key", "Specular Trace Gain Key")
        self.Bind(wx.EVT_MENU, self.on_spec_gain_key, gainThickKey)


        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_view, "&View")
        self.menubar.Append(menu_generate, "&Generate")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)
#v1.2

    def on_press(self,event):
        mylonpt = round(event.xdata,2)
        mylatpt = round(event.ydata,2)
        mylon = round(event.xdata)
        mylat = round(event.ydata)
        press= (mylon,mylat)
        
        if press in self.bullseyes.keys():
            self.bullseyes[press][0].remove() 
            self.bullseyes[press][1].remove() 
            del self.bullseyes[press]  
        else:
            pt = plt.scatter(mylonpt, mylatpt,marker='o',s=9)
            ax = plt.gca()
            annot = ax.annotate('('+str(mylonpt)+','+str(mylatpt)+')',xy=(mylonpt,mylatpt),fontsize=9 )
            self.bullseyes[press] = (pt,annot)
        

    def initialize_plot_map(self):
        ### MAP ###
        screenWidth, screenHeight = wx.GetDisplaySize()
        #print 'W/H',screenWidth,screenHeight
        if screenWidth < 1920 and screenHeight < 1080:
            fig = plt.gcf()
            dpi = fig.dpi
            screenWidthIn = screenWidth/dpi
            screenHeightIn = screenHeight/dpi
            self.height_fig_with_map = 14.0*screenWidthIn/19
            self.width_fig_with_map = 8.0*screenHeightIn/12
        else:
#v1.2
# DPI and figure size was set on development platform kiosk computer
            dpi = 80.0
            self.height_fig_with_map = 14
            self.width_fig_with_map = 8
        background_color = (0/555,76./255,153/255.)
#v1.2
#        self.fig_with_map = plt.figure(figsize=(self.height_fig_with_map, self.width_fig_with_map))
        self.fig_with_map = plt.figure(figsize=(self.height_fig_with_map, self.width_fig_with_map),dpi = dpi)
        self.fig_with_map.set_facecolor(background_color)
        
        
        # Axes
        self.gs_with_map = gridspec.GridSpec(1, 1)
        self.gs_with_map.update(left=0.04, right=0.99, top = 0.98,bottom = 0.04,hspace=0.08,wspace = 0.06)
        self.ax_map = self.fig_with_map.add_subplot(self.gs_with_map[0, 0])

        #This is where we set map params for different zoom settings
        self.zoom_names = ['Globe', 'NWAtlantic', 'SWAtlantic', 'NEPacific', 'SEPacific', 'All Atlantic', 'Asia Pacific', 'Gulf of Mexico', 'Indian Ocean']
        self.min_lon = [-180, -90, -80, -180, -180, -90, 60, -110, 10] #0
        self.max_lon = [180, -10, -0, -100, -100, 30, 180, -30, 140] #60
        self.min_lat = [-90, 0, -50, 0, -50, -40, -40, 0, -70] #-10
        self.max_lat = [90, 50, 0, 50, 0, 40, 40, 50, 40]#40
        self.step_lon = [45,2, 2, 2, 2, 2, 2, 2, 5]#2
        self.step_lat = [30, 2, 2, 2, 2, 2, 2, 2, 5]#2

        #Read in custom zoom levels from ../input_sift/user_zoom.txt
        cNames, cMin_lons, cMax_lons, cMin_lats, cMax_lats, cStep_lons, cStep_lats = readZoom(self.userZoom_filepath)

        #Append info to appropriate list above
        self.zoom_names += cNames
        self.min_lon += cMin_lons
        self.max_lon += cMax_lons
        self.min_lat += cMin_lats
        self.max_lat += cMax_lats
        self.step_lon += cStep_lons
        self.step_lat += cStep_lats

        self.build_plot_map(self.mapRef)

    def build_plot_map(self, mapRef):
        #Del old keys to remove old map lines; must use del here not clear
        for key in self.meridians.keys():
            del self.meridians[key]
        for key in self.parallels.keys():
            del self.parallels[key]

        self.array_lon = [ str(ss) for ss in np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom],self.step_lon[self.zoom]) ]
        self.array_lat = [ str(ss) for ss in np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom],self.step_lat[self.zoom]) ]
        self.ax_map.xaxis.set_major_locator(FixedLocator(np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom]+1, self.step_lon[self.zoom])))
        self.ax_map.yaxis.set_major_locator(FixedLocator(np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom]+1, self.step_lat[self.zoom])))
        self.ax_map.set_xticklabels(self.array_lon, color = 'w', rotation=45)
        self.ax_map.set_yticklabels(self.array_lat, color = 'w')

        self.plotMap = Basemap( projection       = 'cyl',
                 llcrnrlon        = self.min_lon[self.zoom] , #Lower Left  CoRNeR Longitude
                 urcrnrlon        = self.max_lon[self.zoom]  , #Upper Right CoRNeR Longitude
                 llcrnrlat        = self.min_lat[self.zoom]  , #Lower Left  CoRNeR Latitude
                 urcrnrlat        = self.max_lat[self.zoom],   #Upper Right CoRNeR Latitude
                 resolution       = 'l'  ,
                 suppress_ticks   = False,
                 ax = self.ax_map,
                 )
        mapRef = self.plotMap.drawcoastlines(linewidth=0.7, color='blue')
        self.meridians = self.plotMap.drawmeridians(np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom],self.step_lon[self.zoom]))
        self.parallels = self.plotMap.drawparallels(np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom],self.step_lat[self.zoom]))

    def create_main_panel(self):
        self.panel = wx.lib.scrolledpanel.ScrolledPanel(self)
        self.panel.SetupScrolling()
        self.canvas = FigCanvas(self.panel, wx.ID_ANY, self.fig_with_map)

        #IMPORTANT NOTE - SATELLITE NUMBERING
        #Due to the nature of the CYGNSS TLEs
        # - sc 1 is FM05                                                                                                       
        # - sc 2 is FM04                                                                                                       
        # - sc 3 is FM02                                                                                                       
        # - sc 4 is FM01                                                                                                      
        # - sc 5 is FM08                                                                                                       
        # - sc 6 is FM06                                                                                                       
        # - sc 7 is FM07                                                                                                       
        # - sc 8 is FM03
        #This mapping is carried throughout the creation of the main panel
        #However, checkboxes/static text objects are number in terms of sc 1-8
        #This is necessary for SIFT to automatically get the index of a sc from
        #a selected checkbox and map it to the correct satellite object (which are
        #indexed 0-7)

        #Create infrastructure for staticbox outline effect (this must be done first)
        self.paramSb = wx.StaticBox(self.panel, -1, 'Parameters: ')
        self.graphicSb = wx.StaticBox(self.panel, -1, 'SIFT: ')
        self.readoutSb = wx.StaticBox(self.panel, -1, 'Output: ')

        #Create buttons
        self.ppBtn = wx.Button(self.panel, -1, "Play")
        self.Bind(wx.EVT_BUTTON, self.on_pause_button, self.ppBtn)
        self.revBtn = wx.Button(self.panel, -1, "Play Reverse")
        self.Bind(wx.EVT_BUTTON, self.on_rev_button, self.revBtn)

        self.resetBtn = wx.Button(self.panel, -1, "Reset")
        self.Bind(wx.EVT_BUTTON, self.on_reset_button, self.resetBtn)

        self.resetIntervalBtn = wx.Button(self.panel, -1, "Reset Interval")
        self.Bind(wx.EVT_BUTTON, self.on_reset_interval_button, self.resetIntervalBtn)

        self.startValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_startVal_set, self.startValBtn)

        self.endValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_endVal_set, self.endValBtn)

        self.curValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_curVal_set, self.curValBtn)

        self.curPThreeBtn = wx.Button(self.panel, -1, "Set: Current + 3h", name="curPThree", size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.curPThreeBtn)
        self.curPTwelveBtn = wx.Button(self.panel, -1, "Set: Current + 12h", name="curPTwelve", size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.curPTwelveBtn)
        self.realPThreeBtn = wx.Button(self.panel, -1, "Set: Real-Time + 3h", name="realPThree", size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.realPThreeBtn)
        self.realPTwelveBtn = wx.Button(self.panel, -1, "Set: Real-Time + 12h", name="realPTwelve", size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.realPTwelveBtn)

        self.screenshotBtn = wx.Button(self.panel, -1, "Take Screenshot", name="screenshot", size=(120,25))
        self.Bind(wx.EVT_BUTTON, self.on_screenshotBtn, self.screenshotBtn)

       

        #Create static text labels
        self.satLabel = wx.StaticText(self.panel, label="Show/Hide Ground Tracks")
        self.gsLabel = wx.StaticText(self.panel, label="Show/Hide Ground Stations")
        self.speedLabel = wx.StaticText(self.panel, label="Animation Speed")
        self.startLabel = wx.StaticText(self.panel, label="Visualization Start:")
        self.endLabel = wx.StaticText(self.panel, label="Visualization End: ")
        self.jumpLabel = wx.StaticText(self.panel, label="Jump to:                ") #Extra spaces to align textboxes
        self.gsinteractionLabel = wx.StaticText(self.panel, label="Ground Station Interactions ")
        self.stormColorKey = wx.StaticText(self.panel, label="Storm Color Key")

        #For output section
        self.satCol = wx.StaticText(self.panel, label="Satellite")
        self.nextGSCol = wx.StaticText(self.panel, label="Next Ground Station")
        self.nextPassCol = wx.StaticText(self.panel, label="Time until next pass")
        self.durationCol = wx.StaticText(self.panel, label="Duration")

        self.cyg1lab = wx.StaticText(self.panel, label="CYGFM05")
        self.cyg2lab = wx.StaticText(self.panel, label="CYGFM04")
        self.cyg3lab = wx.StaticText(self.panel, label="CYGFM02")
        self.cyg4lab = wx.StaticText(self.panel, label="CYGFM01")
        self.cyg5lab = wx.StaticText(self.panel, label="CYGFM08")
        self.cyg6lab = wx.StaticText(self.panel, label="CYGFM06")
        self.cyg7lab = wx.StaticText(self.panel, label="CYGFM07")
        self.cyg8lab = wx.StaticText(self.panel, label="CYGFM03")
        self.cyg1lab.SetForegroundColour(wx.NamedColour('lawn green')) 
        self.cyg2lab.SetForegroundColour(wx.NamedColour('blue')) 
        self.cyg3lab.SetForegroundColour(wx.NamedColour('purple')) 
        self.cyg4lab.SetForegroundColour(wx.NamedColour('medium orchid')) 
        self.cyg5lab.SetForegroundColour(wx.NamedColour('dodger blue')) 
        self.cyg6lab.SetForegroundColour(wx.NamedColour('steel blue')) 
        self.cyg7lab.SetForegroundColour(wx.NamedColour('sea green')) 
        self.cyg8lab.SetForegroundColour(wx.NamedColour('lime green')) 

        
        self.cyg1NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8NextTC = wx.StaticText(self.panel, size=(100, 20))
      

        self.cygNextTCs = [self.cyg1NextTC, self.cyg2NextTC, self.cyg3NextTC, self.cyg4NextTC, self.cyg5NextTC, self.cyg6NextTC, self.cyg7NextTC, self.cyg8NextTC]

        self.cyg1TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8TimeTC = wx.StaticText(self.panel, size=(100, 20))

        self.cygTimeTCs = [self.cyg1TimeTC, self.cyg2TimeTC, self.cyg3TimeTC, self.cyg4TimeTC, self.cyg5TimeTC, self.cyg6TimeTC, self.cyg7TimeTC, self.cyg8TimeTC]

        self.cyg1DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8DurTC = wx.StaticText(self.panel, size=(100, 20))

        self.cygDurTCs = [self.cyg1DurTC, self.cyg2DurTC, self.cyg3DurTC, self.cyg4DurTC, self.cyg5DurTC, self.cyg6DurTC, self.cyg7DurTC, self.cyg8DurTC]

        #For Storm Color Key
        self.HUlab = wx.StaticText(self.panel, label="Hurricane")
        self.TSlab = wx.StaticText(self.panel, label="Tropical Storm")
        self.TDlab = wx.StaticText(self.panel, label="Tropical Depression")
        self.HUlab.SetForegroundColour(wx.NamedColour(self.satManager.colorByIntensity['HU'])) #red
        self.TSlab.SetForegroundColour(wx.NamedColour(self.satManager.colorByIntensity['TS'])) #orange
        self.TDlab.SetForegroundColour(wx.NamedColour(self.satManager.colorByIntensity['TD'])) #dark turquoise

        #Call to set up readout correctly
        self.updateGSNextAndDur()
        self.updateGSTime()

        #For location information
        self.locationHeader = wx.StaticText(self.panel, label="Instantaneous Satellite Location Information")
        self.satColLab = wx.StaticText(self.panel, label="Satellite")
        self.headingColLab = wx.StaticText(self.panel, label="Heading (deg. from N)")
        self.lonColLab = wx.StaticText(self.panel, label="Longitude")
        self.latColLab = wx.StaticText(self.panel, label="Latitude")
        self.cyg1InfoLab = wx.StaticText(self.panel, label="CYGFM05")
        self.cyg2InfoLab = wx.StaticText(self.panel, label="CYGFM04")
        self.cyg3InfoLab = wx.StaticText(self.panel, label="CYGFM02")
        self.cyg4InfoLab = wx.StaticText(self.panel, label="CYGFM01")
        self.cyg5InfoLab = wx.StaticText(self.panel, label="CYGFM08")
        self.cyg6InfoLab = wx.StaticText(self.panel, label="CYGFM06")
        self.cyg7InfoLab = wx.StaticText(self.panel, label="CYGFM07")
        self.cyg8InfoLab = wx.StaticText(self.panel, label="CYGFM03")
        self.cyg1InfoLab.SetForegroundColour(wx.NamedColour('lawn green'))
        self.cyg2InfoLab.SetForegroundColour(wx.NamedColour('blue')) 
        self.cyg3InfoLab.SetForegroundColour(wx.NamedColour('purple')) 
        self.cyg4InfoLab.SetForegroundColour(wx.NamedColour('medium orchid')) 
        self.cyg5InfoLab.SetForegroundColour(wx.NamedColour('dodger blue')) 
        self.cyg6InfoLab.SetForegroundColour(wx.NamedColour('steel blue')) 
        self.cyg7InfoLab.SetForegroundColour(wx.NamedColour('sea green')) 
        self.cyg8InfoLab.SetForegroundColour(wx.NamedColour('lime green')) 

        self.cyg1HeadTC = wx.StaticText(self.panel)
        self.cyg2HeadTC = wx.StaticText(self.panel)
        self.cyg3HeadTC = wx.StaticText(self.panel)
        self.cyg4HeadTC = wx.StaticText(self.panel)
        self.cyg5HeadTC = wx.StaticText(self.panel)
        self.cyg6HeadTC = wx.StaticText(self.panel)
        self.cyg7HeadTC = wx.StaticText(self.panel)
        self.cyg8HeadTC = wx.StaticText(self.panel)

        self.cygHeadTCs = [self.cyg1HeadTC, self.cyg2HeadTC, self.cyg3HeadTC, self.cyg4HeadTC, self.cyg5HeadTC, self.cyg6HeadTC, self.cyg7HeadTC, self.cyg8HeadTC]

        self.cyg1LonTC = wx.StaticText(self.panel)
        self.cyg2LonTC = wx.StaticText(self.panel)
        self.cyg3LonTC = wx.StaticText(self.panel)
        self.cyg4LonTC = wx.StaticText(self.panel)
        self.cyg5LonTC = wx.StaticText(self.panel)
        self.cyg6LonTC = wx.StaticText(self.panel)
        self.cyg7LonTC = wx.StaticText(self.panel)
        self.cyg8LonTC = wx.StaticText(self.panel)

        self.cygLonTCs = [self.cyg1LonTC, self.cyg2LonTC, self.cyg3LonTC, self.cyg4LonTC, self.cyg5LonTC, self.cyg6LonTC, self.cyg7LonTC, self.cyg8LonTC]

        self.cyg1LatTC = wx.StaticText(self.panel)
        self.cyg2LatTC = wx.StaticText(self.panel)
        self.cyg3LatTC = wx.StaticText(self.panel)
        self.cyg4LatTC = wx.StaticText(self.panel)
        self.cyg5LatTC = wx.StaticText(self.panel)
        self.cyg6LatTC = wx.StaticText(self.panel)
        self.cyg7LatTC = wx.StaticText(self.panel)
        self.cyg8LatTC = wx.StaticText(self.panel)

        self.cygLatTCs = [self.cyg1LatTC, self.cyg2LatTC, self.cyg3LatTC, self.cyg4LatTC, self.cyg5LatTC, self.cyg6LatTC, self.cyg7LatTC, self.cyg8LatTC]

        #Set initial headings
        self.updateLocAndHead()

        self.vizParamLab = wx.StaticText(self.panel, size=(150,20), label="Visualization Parameters")
        self.startParamLab = wx.StaticText(self.panel, size=(120,20), label="Current Analysis Start: ")
        self.endParamLab = wx.StaticText(self.panel, size=(120,20), label="Current Analysis End: ")
        self.maxStartParamLab = wx.StaticText(self.panel, size=(100,20), label="Max Analysis Start: ")
        self.maxEndParamLab = wx.StaticText(self.panel, size=(100,20), label = "Max Analysis End: ")

        self.startParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.start_visu)
        self.endParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.stop_visu)
        self.maxStartParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.firstTime)
        self.maxEndParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.lastTime)

        #Create checkboxes
        self.cyg1cb = wx.CheckBox(self.panel, id=1, label='CYGFM05')
        self.cyg1cb.SetForegroundColour(wx.NamedColour('lawn green')) 
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg1cb)

        self.cyg2cb = wx.CheckBox(self.panel, id=2, label='CYGFM04')
        self.cyg2cb.SetForegroundColour(wx.NamedColour('blue'))
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg2cb)

        self.cyg3cb = wx.CheckBox(self.panel, id=3, label='CYGFM02')
        self.cyg3cb.SetForegroundColour(wx.NamedColour('purple4')) 
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg3cb)

        self.cyg4cb = wx.CheckBox(self.panel, id=4, label='CYGFM01',)
        self.cyg4cb.SetForegroundColour(wx.NamedColour('medium orchid')) 
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg4cb)

        self.cyg5cb = wx.CheckBox(self.panel, id=5, label='CYGFM08')
        self.cyg5cb.SetForegroundColour(wx.NamedColour('dodger blue')) 
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg5cb)

        self.cyg6cb = wx.CheckBox(self.panel, id=6, label='CYGFM06')
        self.cyg6cb.SetForegroundColour(wx.NamedColour('steel blue'))
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg6cb)

        self.cyg7cb = wx.CheckBox(self.panel, id=7, label='CYGFM07')
        self.cyg7cb.SetForegroundColour(wx.NamedColour('sea green'))
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg7cb)

        self.cyg8cb = wx.CheckBox(self.panel, id=8, label='CYGFM03')
        self.cyg8cb.SetForegroundColour(wx.NamedColour('lime green'))
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg8cb)

        self.cyg1specCB = wx.CheckBox(self.panel, id=9, label='Sp 5')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg1specCB)
        self.cyg2specCB = wx.CheckBox(self.panel, id=10, label='Sp 4')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg2specCB)
        self.cyg3specCB = wx.CheckBox(self.panel, id=11, label='Sp 2')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg3specCB)
        self.cyg4specCB = wx.CheckBox(self.panel, id=12, label='Sp 1')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg4specCB)
        self.cyg5specCB = wx.CheckBox(self.panel, id=13, label='Sp 8')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg5specCB)
        self.cyg6specCB = wx.CheckBox(self.panel, id=14, label='Sp 6')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg6specCB)
        self.cyg7specCB = wx.CheckBox(self.panel, id=15, label='Sp 7')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg7specCB)
        self.cyg8specCB = wx.CheckBox(self.panel, id=16, label='Sp 3')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg8specCB)

        # self.specGainCB = wx.CheckBox(self.panel, label='Display Sp Gain (slow)')
        # self.Bind(wx.EVT_CHECKBOX, self.on_cb_specGain, self.specGainCB)

        self.cyg1gainCB = wx.CheckBox(self.panel, id=17, label='G 5')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg1gainCB)
        self.cyg2gainCB = wx.CheckBox(self.panel, id=18, label='G 4')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg2gainCB)
        self.cyg3gainCB = wx.CheckBox(self.panel, id=19, label='G 2')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg3gainCB)
        self.cyg4gainCB = wx.CheckBox(self.panel, id=20, label='G 1')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg4gainCB)
        self.cyg5gainCB = wx.CheckBox(self.panel, id=21, label='G 8')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg5gainCB)
        self.cyg6gainCB = wx.CheckBox(self.panel, id=22, label='G 6')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg6gainCB)
        self.cyg7gainCB = wx.CheckBox(self.panel, id=23, label='G 7')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg7gainCB)
        self.cyg8gainCB = wx.CheckBox(self.panel, id=24, label='G 3')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_gain, self.cyg8gainCB)

        self.gs1cb = wx.CheckBox(self.panel, label='Hawaii')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsHawaii, self.gs1cb)
        self.gs2cb = wx.CheckBox(self.panel, label='Chile')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsSantiago, self.gs2cb)
        self.gs3cb = wx.CheckBox(self.panel, label='Australia')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsPerth, self.gs3cb)

        self.stormcb = wx.CheckBox(self.panel, label='Show/Hide Storms')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_storms, self.stormcb)

        self.cygCheckBoxes = [self.cyg1cb, self.cyg2cb, self.cyg3cb, self.cyg4cb, self.cyg5cb, self.cyg6cb, self.cyg7cb, self.cyg8cb]
        self.gsCheckBoxes = [self.gs1cb, self.gs2cb, self.gs3cb]
        self.specCheckBoxes = [self.cyg1specCB, self.cyg2specCB, self.cyg3specCB, self.cyg4specCB, self.cyg5specCB, self.cyg6specCB, self.cyg7specCB, self.cyg8specCB]

        #Create sliders
        self.speedSld = wx.Slider(self.panel, value=self.satManager.speedFactor, minValue=self.satManager.speedFactor, maxValue=100, style=wx.SL_HORIZONTAL)
        self.speedSld.Bind(wx.EVT_SCROLL, self.on_speed_scroll)

        #Create text entry fields
        self.startTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(100,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_startVal_set, self.startTC)
        self.endTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(100,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_endVal_set, self.endTC)
        self.curTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(100,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_curVal_set, self.curTC)


        #Main Sizers
        self.windowSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.controlsSizer = wx.StaticBoxSizer(self.paramSb, wx.VERTICAL)
        self.outputSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.graphicSizer = wx.StaticBoxSizer(self.graphicSb, wx.VERTICAL)
        self.textoutSizer = wx.StaticBoxSizer(self.readoutSb, wx.VERTICAL)

        #Row sizers for controls section
        self.cRow1 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow2 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow3 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow4 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow5 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow6 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow7 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow8 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow9 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow10 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow11 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow12 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow13 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow14 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow15 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow16 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow17 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow18 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow19 = wx.BoxSizer(wx.HORIZONTAL)

        #Row sizers for output section
        self.oRow0 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow1 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow2 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow3 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow4 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow5 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow6 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow7 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow8 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow9 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow10 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow11 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow12 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow13 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow14 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow15 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow16 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow17 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow18 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow19 = wx.BoxSizer(wx.HORIZONTAL)
        #(Row, Col, vgap, hgap)
        self.satIdSizer = wx.GridSizer(9,1,2,2)
        self.headColSizer = wx.GridSizer(9,1,2,2)
        self.lonColSizer = wx.GridSizer(9,1,2,2)
        self.latColSizer = wx.GridSizer(9,1,2,2)

        #Fill sizers in order from smallest to largest
        #Control rows
        #Note the apparent out of order cb adding yields an in order listing
        #For details, see IMPORTANT NOTE - SATELLITE NUMBERING (use ctrl+f)
        self.cRow1.Add(self.cyg4cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg4specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg1cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg1specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow2.Add(self.cyg3cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg3specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg6cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg6specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow3.Add(self.cyg8cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg8specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg7cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg7specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow4.Add(self.cyg2cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg2specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg5cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg5specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow17.Add(self.cyg4gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow17.Add(self.cyg3gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow17.Add(self.cyg8gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow17.Add(self.cyg2gainCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow19.Add(self.cyg1gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow19.Add(self.cyg6gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow19.Add(self.cyg7gainCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow19.Add(self.cyg5gainCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow5.Add(self.gs1cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow6.Add(self.gs2cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow7.Add(self.gs3cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow8.Add(self.stormcb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow9.Add(self.speedLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow9.Add(self.speedSld, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow10.Add(self.startLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow10.Add(self.startTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow10.Add(self.startValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow11.Add(self.endLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow11.Add(self.endTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow11.Add(self.endValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow12.Add(self.jumpLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow12.Add(self.curTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow12.Add(self.curValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow13.Add(self.ppBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow13.Add(self.revBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow16.Add(self.resetBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow16.Add(self.resetIntervalBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow14.Add(self.curPThreeBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow14.Add(self.curPTwelveBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow15.Add(self.realPThreeBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow15.Add(self.realPTwelveBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow18.Add(self.screenshotBtn, 0, wx.ALL | wx.EXPAND, border = 4)


        #Output rows
        self.oRow0.Add(self.satCol, 0, wx.ALL | wx.EXPAND, border= 4)
        self.oRow0.Add(self.nextGSCol, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow0.Add(self.nextPassCol, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow0.Add(self.durationCol, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow1.Add(self.cyg4lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg4NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg4TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg4DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow2.Add(self.cyg3lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg3NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg3TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg3DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow3.Add(self.cyg8lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg8NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg8TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg8DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow4.Add(self.cyg2lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg2NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg2TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg2DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow5.Add(self.cyg1lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg1NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg1TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg1DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow6.Add(self.cyg6lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow7.Add(self.cyg7lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow8.Add(self.cyg5lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg5NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg5TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg5DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow9.Add(self.vizParamLab, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow10.Add(self.startParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow10.Add(self.startParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow11.Add(self.endParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow11.Add(self.endParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow12.Add(self.maxStartParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow12.Add(self.maxStartParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow13.Add(self.maxEndParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow13.Add(self.maxEndParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow14.Add(self.stormColorKey, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow15.Add(self.HUlab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow16.Add(self.TSlab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow17.Add(self.TDlab, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow18.Add(self.locationHeader, 0, wx.ALL | wx.EXPAND, border = 4)

        self.satIdSizer.Add(self.satColLab, 0, wx.ALL | wx.EXPAND, border=4)
        self.satIdSizer.Add(self.cyg4InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg3InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg8InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg2InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg1InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg6InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg7InfoLab, 0, wx.ALL, border=4)
        self.satIdSizer.Add(self.cyg5InfoLab, 0, wx.ALL, border=4)

        self.headColSizer.Add(self.headingColLab, 0, wx.ALL | wx.EXPAND, border=4)
        self.headColSizer.Add(self.cyg4HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg3HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg8HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg2HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg1HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg6HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg7HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.headColSizer.Add(self.cyg5HeadTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)

        self.lonColSizer.Add(self.lonColLab, 0, wx.ALL | wx.EXPAND, border=4)
        self.lonColSizer.Add(self.cyg4LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg3LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg8LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg2LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg1LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg6LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg7LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.lonColSizer.Add(self.cyg5LonTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        
        self.latColSizer.Add(self.latColLab, 0, wx.ALL | wx.EXPAND, border=4)
        self.latColSizer.Add(self.cyg4LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)   
        self.latColSizer.Add(self.cyg3LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg8LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg2LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg1LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg6LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg7LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)
        self.latColSizer.Add(self.cyg5LatTC, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=4)

        self.oRow19.Add(self.satIdSizer, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow19.Add(self.headColSizer, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow19.Add(self.lonColSizer, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow19.Add(self.latColSizer, 0, wx.ALL | wx.EXPAND, border=4)

        self.controlsSizer.Add(self.satLabel, 0, wx.ALL, border = 10)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow1, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow2, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow3, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow4, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow17, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow19, 0, wx.ALL | wx.EXPAND, border = 2)

        self.controlsSizer.Add(self.gsLabel, 0, wx.ALL, border = 10)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow5, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow6, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow7, 0, wx.ALL | wx.EXPAND, border = 2)

        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow8, 0, wx.ALL | wx.EXPAND, border = 2)

        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow9, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow10, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow11, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow12, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow13, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow16, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow14, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow15, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow18, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)


        self.textoutSizer.Add(self.gsinteractionLabel, 0, wx.ALL, border=5)
        self.textoutSizer.Add(self.oRow0, 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow1, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow2, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow3, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow4, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow5, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow6, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow7, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow8, 0, wx.ALL | wx.EXPAND, border=2)

        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow18, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow19, 0, wx.ALL | wx.EXPAND, border=2)

        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow9, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow10, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow11, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow12, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow13, 0, wx.ALL | wx.EXPAND, border=2)

        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow14, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow15, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow16, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow17, 0, wx.ALL | wx.EXPAND, border=2)

        self.graphicSizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND, border = 15)
        #self.graphicSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        #self.graphicSizer.Add(self.outputLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        #self.graphicSizer.Add(self.gRow0, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)

        # self.outputSizer.Add(self.graphicSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        # self.outputSizer.Add(self.textoutSizer, 0, wx.ALL | wx.EXPAND, border = 5)

        self.windowSizer.Add(self.controlsSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        # self.windowSizer.Add(self.outputSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        self.windowSizer.Add(self.graphicSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        self.windowSizer.Add(self.textoutSizer, 0, wx.ALL | wx.EXPAND, border = 5)

#        self.controlsSizer.Fit(self)
#        self.outputSizer.Fit()
#        self.graphicSizer.Fit(self)
#        self.textoutSizer.Fit(self)

#v1.2
        # not commented self.panel.SetSizer(self.windowSizer)
        # panelSizer = wx.BoxSizer(wx.HORIZONTAL)
        # panelSizer.Add(self.panel)
        # self.SetSizerAndFit(panelSizer)
        # not commented self.windowSizer.Fit(self)

        self.panel.SetSizer(self.windowSizer)
        # panelSizer = wx.BoxSizer(wx.HORIZONTAL)
        # panelSizer.Add(self.panel)
        # self.SetSizerAndFit(self.windowSizer)
        self.windowSizer.Fit(self)

    #Helper functions ###################################

    def draw_plot(self):
        self.satManager.advSim(self.plotMap)
        self.updateGSNextAndDur()
        self.updateGSTime()
        self.updateLocAndHead()

        self.canvas.draw()

    def on_redraw_timer(self, event):
        if not self.paused:
            self.satManager.incrementIndex()
        elif not self.revPaused:
            self.satManager.decrementIndex()

        self.draw_plot()


    #GUI Control Functionality ###########################    

    def on_pause_button(self, event):
        self.paused = not self.paused

        label = "Play" if self.paused else "Pause"
        self.ppBtn.SetLabel(label)

        if self.paused:
            self.revBtn.Enable()
            self.satManager.animDirection = 0
        else:
            self.revBtn.Disable()
            self.satManager.animDirection = 1

    def on_rev_button(self, event):
        self.revPaused = not self.revPaused

        label = "Play Reversed" if self.revPaused else "Pause"
        self.revBtn.SetLabel(label)

        if self.revPaused:
            self.ppBtn.Enable()
            self.satManager.animDirection = 0
        else:
            self.ppBtn.Disable()
            self.satManager.animDirection = -1

    def on_reset_button(self, event):
        #We choose to automatically pause when resetting
        self.paused = True
        self.ppBtn.SetLabel("Play")

        self.satManager.curSample = self.satManager.when_start_visu

        #Perform updates to GS indexing to ensure readout is correct
        self.satManager.findGSinteractionIdx()
        self.updateGSNextAndDur()
        self.updateGSTime()

        self.satManager.removeStorms()
        self.satManager.updateCurStorm()
        self.satManager.plotStorms(self.plotMap)

    def on_reset_interval_button(self, event):
        self.satManager.start_visu = self.satManager.firstTime
        self.satManager.stop_visu = self.satManager.lastTime
        #These values are indices
        self.satManager.when_start_visu = self.satManager.time_sat.index(self.satManager.start_visu)
        self.satManager.when_stop_visu = self.satManager.time_sat.index(self.satManager.stop_visu)

        #Update visualization parameter readout
        self.startParam.SetLabel(self.satManager.firstTime)
        self.endParam.SetLabel(self.satManager.lastTime)

        logging.info("New Analysis Interval: %s -> %s", self.satManager.time_sat[self.satManager.when_start_visu], self.satManager.time_sat[self.satManager.when_stop_visu])

    def on_hotBtn(self, event):
        caller = event.GetEventObject().GetName()
        
        if caller == "curPThree":
            newStartTime = datetime.strptime(self.satManager.time_sat[self.satManager.curSample],'%Y-%m-%dT%H:%M:%S')
            newEndTime = newStartTime + timedelta(hours=3)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

        elif caller == "curPTwelve":
            newStartTime = datetime.strptime(self.satManager.time_sat[self.satManager.curSample],'%Y-%m-%dT%H:%M:%S')
            newEndTime = newStartTime + timedelta(hours=12)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

        elif caller == "realPThree":
            newStartTime = datetime.utcnow()
            newEndTime = newStartTime + timedelta(hours=3)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))


        elif caller == "realPTwelve":
            newStartTime = datetime.utcnow()
            newEndTime = newStartTime + timedelta(hours=12)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

    def on_screenshotBtn(self, event):
        screen = wx.ScreenDC()
        size = screen.GetSize()
        bmp = wx.EmptyBitmap(size[0], size[1])
        mem = wx.MemoryDC(bmp)
        mem.Blit(0, 0, size[0], size[1], screen, 0, 0)
        del mem  # Release bitmap

        if not os.path.isdir(self.satManager.screenshot_filepath):
            os.makedirs(self.satManager.screenshot_filepath)

        timestamp = time.strftime('%Y-%m-%dT%H-%M-%S')
        bmp.SaveFile(self.satManager.screenshot_filepath + '/screenshot_' + timestamp +'.png', wx.BITMAP_TYPE_PNG)
        logging.info("Screenshot %s.png captured", timestamp)


    def time_set_helper(self, newStartTime, newEndTime):
        try:
            indexStartTime = self.satManager.time_sat.index(newStartTime)
            indexEndTime = self.satManager.time_sat.index(newEndTime)
            #If there is no exception in the above statement, set values
            #Note: start_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_start_visu is index corresponding to that timestamp
            self.satManager.start_visu = newStartTime
            self.satManager.when_start_visu = indexStartTime
            self.satManager.curSample = indexStartTime

            self.satManager.stop_visu = newEndTime
            self.satManager.when_stop_visu = indexEndTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate GS interaction indices
            self.satManager.findGSinteractionIdx()

            #Update visualization param readout
            self.startParam.SetLabel(newStartTime)
            self.endParam.SetLabel(newEndTime)
            
        except ValueError:
            logging.error("OPERATION FAILED: %s -> %s is outside allowable range; Analysis start time has not changed", newStartTime, newEndTime)
            return

        logging.info("New Analysis Interval: %s -> %s", self.satManager.time_sat[self.satManager.when_start_visu], self.satManager.time_sat[self.satManager.when_stop_visu])


    def on_cb_cyg(self, event):
        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        senderLabel = sender.GetLabel()

        #Sat cb are numbered 1 - 8, subtract 1 for zero-indexing
        cbId = int(sender.GetId()) - 1

        if(isChecked):
            self.satManager.plotSatTrace(self.plotMap, cbId)
        else:
            self.satManager.removeSatTrace(cbId)

    def on_cb_cyg_spec(self, event):
        sender = event.GetEventObject()
        isChecked=sender.GetValue()

        senderLabel = sender.GetLabel()

        #Spec cb are numbered 8 to 16, so we subtract offset of 8 and 1 for zero-indexing
        cbId = int(sender.GetId()) - 8 - 1

        if(isChecked):
            self.satManager.plotSpecTrace(self.plotMap, cbId)
        else:
            self.satManager.removeSpecTrace(cbId)

    def on_cb_cyg_gain(self, event):
        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        senderLabel = sender.GetLabel()

        #Gain cb are numbered 17 to 24, so subtract offset of 16 and 1 for zero-indexing
        cbId = int(sender.GetId()) - 16 - 1

        if(isChecked):
            self.satManager.plotSpecGainTrace(self.plotMap, cbId)
        else:
            self.satManager.removeSpecGainTrace(cbId)


    def on_cb_gsHawaii(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 0)
        else:
            self.satManager.removeGroundStation(0)

#        self.draw_plot()

    def on_cb_gsSantiago(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 1)
        else:
            self.satManager.removeGroundStation(1)

    def on_cb_gsPerth(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 2)
        else:
            self.satManager.removeGroundStation(2)


    def on_cb_storms(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotTrajectory(self.plotMap)
            self.satManager.plotStorms(self.plotMap)
        else:
            self.satManager.removeTrajectory()
            self.satManager.removeStorms()

    def on_speed_scroll(self, event):

        sender = event.GetEventObject()
        newFactor = sender.GetValue()

        self.satManager.speedFactor = newFactor

    def on_startVal_set(self, event):
        newTime = self.startTC.GetValue()

        #Check if input time is a valid timestamp
        try:
            indexTime = self.satManager.time_sat.index(newTime)
            #If there is no exception in the above statement, set values
            #Note: start_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_start_visu is index corresponding to that timestamp
            self.satManager.start_visu = newTime
            self.satManager.when_start_visu = indexTime
            self.satManager.curSample = indexTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate GS interaction indices
            self.satManager.findGSinteractionIdx()

            #Recalculate spec/storm interaction
            #self.satManager.findSpecStormPasses()

            #Update visualization param readout
            self.startParam.SetLabel(newTime)

            logging.info("New Analysis Interval: %s -> %s", self.satManager.time_sat[self.satManager.when_start_visu], self.satManager.time_sat[self.satManager.when_stop_visu])
    
        except ValueError:
            logging.error("OPERATION FAILED: Start time entered is outside allowable range; Analysis start time has not changed.")

        
            
    def on_endVal_set(self, event):
        newTime = self.endTC.GetValue()

        #Check if input time is a valid timestamp
        try:
            indexTime = self.satManager.time_sat.index(newTime)
            #If there is no exception in the above statement, set values
            #Note: stop_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_stop_visu is index corresponding to that timestamp
            self.satManager.stop_visu = newTime
            self.satManager.when_stop_visu = indexTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate spec/storm interaction
            #self.satManager.findSpecStormPasses()
            #print self.satManager.specPassDist

            self.endParam.SetLabel(newTime)

            logging.info("New Analysis Interval: %s -> %s", self.satManager.time_sat[self.satManager.when_start_visu], self.satManager.time_sat[self.satManager.when_stop_visu])
    
        except ValueError:
            logging.error("OPERATION FAILED: End time entered is outside allowable range; Analysis end time has not changed.")

            #Update visualization param readout
            self.endParam.SetLabel(newTime)

        

    def on_curVal_set(self, event):
        newTime = self.curTC.GetValue()

        #First check to make sure input is a valid value at all
        try:
            curTimeIdx = self.satManager.time_sat.index(newTime)
            
        except ValueError:
            logging.error("OPERATION FAILED: Invalid timestamp entered; Please enter a timestamp between %s and %s.", self.satManager.start_visu, self.satManager.stop_visu)
            return

        #Now test if input is between CURRENT start and stop times
        if curTimeIdx <= self.satManager.when_stop_visu and curTimeIdx >= self.satManager.when_start_visu:
            self.satManager.curSample = curTimeIdx
        else:
            logging.error("OPERATION FAILED: Invalid timestamp entered; Please enter a timestamp between %s and %s.", self.satManager.start_visu, self.satManager.stop_visu)


    def on_exit(self, event):
        self.Destroy()
        self.redraw_timer.Stop()
        sys.exit(0)


    def on_zoom(self, event):
        #Get info about what menu item was clicked
        menubar = event.GetEventObject().GetMenuBar()
        id_selected = event.GetId()
        menuItem = menubar.FindItemById(id_selected).GetLabel()

        #Use info to index into zoom array
        self.zoom = self.zoom_names.index(menuItem)
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_overpass_report(self, event):
        areaSelectWindow = OverpassReportFrame(self.satManager)
        areaSelectWindow.Show()

    def on_generate_spec_report(self, event):
        self.satManager.generateSpecReport()

    def on_generate_gs_report(self, event):
        try:
            self.satManager.generateGSreport()
            successMsg = wx.MessageDialog(self, "CYGNSS ground station report generated.", "CYGNSS Ground Station Report", style=wx.OK|wx.CENTRE)
            successMsg.ShowModal()
        except Exception, e:
            logging.error(str(e))
            failureMsg = wx.MessageDialog(self, "CYGNSS ground station report failed.", "CYGNSS Ground Station Report", style=wx.OK|wx.CENTRE)
            failureMsg.ShowModal()


    def on_storm_color_key(self, event):
        colorKeyWindow = ColorKeyFrame(self.satManager.colorByIntensity, self.satManager.colorByIntensity_faded)
        colorKeyWindow.Show()

    def on_spec_gain_key(self, event):
        specGainWindow = GainKeyFrame()
        specGainWindow.Show()


    def update_map(self):
        self.initialize_plot_map()
        #self.satManager = SatelliteManager(self.plotMap, self.max_lon, self.min_lat, self.height_fig_with_map, self.width_fig_with_map, self.ax_map)
        self.canvas = FigCanvas(self.panel, wx.ID_ANY, self.fig_with_map)

    def reset_sat_cb(self):
        for i in range(len(self.cygCheckBoxes)):
            self.cygCheckBoxes[i].SetValue(False)

    def reset_traces(self):
        for idx in range(len(self.cygCheckBoxes)):
            #print self.cygCheckBoxes[idx].IsChecked()
            if self.cygCheckBoxes[idx].IsChecked():
                self.satManager.removeSatTrace(idx)
                self.satManager.plotSatTrace(self.plotMap, idx)

            if self.specCheckBoxes[idx].IsChecked():
                self.satManager.removeSpecTrace(idx)
                self.satManager.plotSpecTrace(self.plotMap, idx)

    #GUI Info Readout Functionality #######################

    def updateGSNextAndDur(self):
        for i in range(self.satManager.nb_satellites):
            upperLimit = len(self.satManager.overGSstart[i])
            lowerLimit = -1
            if(self.satManager.GSinteractionIdx[i] != self.localGSInteractionIdxs[i]):
                self.localGSInteractionIdxs[i] = self.satManager.GSinteractionIdx[i]
                if self.satManager.GSinteractionIdx[i] < upperLimit and self.satManager.GSinteractionIdx[i] > lowerLimit:
                    self.cygNextTCs[i].SetLabel(self.satManager.overGSwhich[i][self.satManager.GSinteractionIdx[i]])
                    self.cygDurTCs[i].SetLabel(str(self.satManager.overGSdur[i][self.satManager.GSinteractionIdx[i]]))
                else: #Next GS is outside analysis interval
                    self.cygNextTCs[i].SetLabel("None")
                    self.cygDurTCs[i].SetLabel("-")
        

    def updateGSTime(self):
        for i in range(self.satManager.nb_satellites):
            upperLimit = len(self.satManager.overGSstart[i])
            lowerLimit = -1
            if self.satManager.GSinteractionIdx[i] < upperLimit and self.satManager.GSinteractionIdx[i] > lowerLimit:
                time = (self.satManager.overGSstart[i][self.satManager.GSinteractionIdx[i]]) - datetime.strptime(self.satManager.time_sat[self.satManager.curSample], '%Y-%m-%dT%H:%M:%S')
                self.cygTimeTCs[i].SetLabel(str(time))
            else: #Next GS is outside analysis interval
                self.cygTimeTCs[i].SetLabel('-')

    def updateLocAndHead(self):
        for i in range(self.satManager.nb_satellites):
            self.cygHeadTCs[i].SetLabel(str(round(self.satManager.heading_sat[i,self.satManager.curSample],0)))
            self.cygLonTCs[i].SetLabel(str(round(self.satManager.lon_sat[i, self.satManager.curSample],2)))
            self.cygLatTCs[i].SetLabel(str(round(self.satManager.lat_sat[i, self.satManager.curSample],2)))

class ColorKeyFrame(wx.Frame):
    title = 'Forecast Color Key'

    def __init__(self, colors, fadedColors):
        wx.Frame.__init__(self, None, -1, self.title)
        self.panel = wx.Panel(self, -1)
        self.SetBackgroundColour("light gray")

        #Map storm type shortcodes to human readable names
        shortcodeMap = {'DB':'Disturbance', 'TD':'Tropical Depression', 'TS':'Tropical Storm', 'TY':'Typhoon', 'ST':'Super Typhoon', 'TC':'Tropical Cyclone',
                                    'HU':'Hurricane', 'SD':'Subtropical Depression', 'SS':'Subtropical Storm', 'EX':'Extratropical Systems', 'PT':'Post Tropical',
                                    'IN':'Inland', 'DS':'Dissipating', 'LO': 'Low', 'WV': 'Tropical Wave', 'ET':'Extrapolated',
                                    'XX':'Unknown'}

        vsizer = wx.BoxSizer(wx.VERTICAL)
        # wx.GridSizer(rows, cols, vgap, hgap)
        gsizer = wx.GridSizer(len(colors), 3, 2, 2)

        for key, value in colors.items():
            name = wx.StaticText(self, wx.ID_ANY, label=shortcodeMap[key])
            gsizer.Add(name, 0, wx.ALL, border=2)
            
            color = wx.Panel(self, wx.ID_ANY)
            color.SetBackgroundColour(wx.NamedColour(value))
            gsizer.Add(color, 0, wx.ALL|wx.EXPAND, border=2)

            fadedColor = wx.Panel(self, wx.ID_ANY)
            fadedColor.SetBackgroundColour(wx.NamedColour(fadedColors[key]))
            gsizer.Add(fadedColor, 0, wx.ALL|wx.EXPAND, border=2)

        vsizer.Add(gsizer, 0, wx.ALL|wx.EXPAND, 10)

        self.SetSizerAndFit(vsizer)
        self.Centre()
        self.Show()

class GainKeyFrame(wx.Frame):
    title = 'Specular Trace Gain Thickness Key'

    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        panel = wx.Panel(self, -1)
        fig = plt.figure(figsize=(3,3))
        fig.patch.set_facecolor('whitesmoke')
        ax = fig.add_subplot(111)

        keySizer = wx.BoxSizer(wx.HORIZONTAL)

        x = [0,1]
        highy = [3,3]
        medy = [2,2]
        lowy = [1,1]

        ax.scatter(x, highy, s=15, c='darkgrey')
        ax.plot(x, highy, linewidth=5, color='darkgrey')
        ax.text(1.2,3, "High Gain (>10 dB)", verticalalignment='center')
        ax.scatter(x, medy, s=5, c='darkgrey')
        ax.plot(x, medy, linewidth=3.5, color='darkgrey')
        ax.text(1.2,2, "Med Gain (>5 dB)", verticalalignment='center')
        ax.scatter(x, lowy, s=0.5, c='darkgrey')
        ax.plot(x, lowy, linewidth=0.9, color='darkgrey')
        ax.text(1.2,1, "Low Gain (<5 dB)", verticalalignment='center')

        ax.set_xlim([-0.1,2])
        ax.set_ylim([0,3.5])
        
        ax.axis('off')

        canvas = FigCanvas(panel, wx.ID_ANY, fig)
        canvas.draw()

        keySizer.Add(canvas, 1, wx.ALL | wx.EXPAND, border = 0)
        panel.SetSizerAndFit(keySizer)

class OverpassReportFrame(wx.Frame):
    title = 'CYGNSS Overpass Report'
    satManager = None

    def __init__(self, inSatManager):
        style = (wx.MINIMIZE_BOX | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX |
                 wx.CLIP_CHILDREN | wx.FRAME_FLOAT_ON_PARENT)
        wx.Frame.__init__(self, wx.GetApp().TopWindow, -1, self.title, style=style, size=(300,300))
        self.Bind(wx.EVT_CLOSE, self.on_close)
        panel = wx.Panel(self, -1)
        self.SetBackgroundColour("light gray")
        # self.MakeModal()

        self.satManager = inSatManager

        #Create sizers
        vsizer = wx.BoxSizer(wx.VERTICAL)
        cbSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        gsizer = wx.GridSizer(6,2,2,2)

        #Create checkbox
        self.cbPoints = wx.CheckBox(panel, label="Include sample locations?")

        #Create Labels
        llLonLab = wx.StaticText(panel, label="Lower Left Lon: ")
        llLatLab = wx.StaticText(panel, label="Lower Left Lat: ")
        urLonLab = wx.StaticText(panel, label="Upper Right Lon: ")
        urLatLab = wx.StaticText(panel, label="Upper Right Lat: ")
        startLab = wx.StaticText(panel, label="Report Start Time: ")
        endLab = wx.StaticText(panel, label="Report End Time: ")

        #Create Text Ctrls
        self.llLonTC = wx.TextCtrl(panel)
        self.llLatTC = wx.TextCtrl(panel)
        self.urLonTC = wx.TextCtrl(panel)
        self.urLatTC = wx.TextCtrl(panel)
        self.startTC = wx.TextCtrl(panel)
        self.endTC = wx.TextCtrl(panel)

        #Create Buttons
        okButton = wx.Button(panel, label='OK')
        self.Bind(wx.EVT_BUTTON, self.on_ok_button, okButton)
        cancelButton = wx.Button(panel, label='Cancel')
        self.Bind(wx.EVT_BUTTON, self.on_close, cancelButton)

        #Assemble
        cbSizer.Add(self.cbPoints, 0, wx.ALL, border=2)

        gsizer.Add(llLonLab, 0, wx.ALL, border=2)
        gsizer.Add(self.llLonTC, 0, wx.ALL, border=2)
        gsizer.Add(llLatLab, 0, wx.ALL, border=2)
        gsizer.Add(self.llLatTC, 0, wx.ALL, border=2)
        gsizer.Add(urLonLab, 0, wx.ALL, border=2)
        gsizer.Add(self.urLonTC, 0, wx.ALL, border=2)
        gsizer.Add(urLatLab, 0, wx.ALL, border=2)
        gsizer.Add(self.urLatTC, 0, wx.ALL, border=2)
        gsizer.Add(startLab, 0, wx.ALL, border=2)
        gsizer.Add(self.startTC, 0, wx.ALL, border=2)
        gsizer.Add(endLab, 0, wx.ALL, border=2)
        gsizer.Add(self.endTC, 0, wx.ALL, border=2)

        buttonSizer.Add(okButton, 0, wx.EXPAND|wx.ALL, border=4)
        buttonSizer.Add(cancelButton, 0, wx.EXPAND|wx.ALL, border=4)

        vsizer.Add(cbSizer, 0, wx.ALL|wx.EXPAND, 10)
        vsizer.Add(gsizer)
        vsizer.Add(buttonSizer)

        panel.SetSizer(vsizer)
        self.Centre()

    def on_ok_button(self, event):
        #Check for values in all boxes and if valid, call report generator
        llLat = self.llLatTC.GetValue()
        llLon = self.llLonTC.GetValue()
        urLat = self.urLatTC.GetValue()
        urLon = self.urLonTC.GetValue()
        try:
            start = datetime.strptime(self.startTC.GetValue(), '%Y-%m-%dT%H:%M:%S')
            end = datetime.strptime(self.endTC.GetValue(), '%Y-%m-%dT%H:%M:%S')
        except ValueError:
            errorMsg = wx.MessageDialog(self, "Error: Invalid datetime format - use YYYY-mm-DDTHH:MM:SS", "Error", style=wx.OK|wx.CENTRE)
            errorMsg.ShowModal()
            return

        #Check times are within analysis interval bounds
        if not all (x >= datetime.strptime(self.satManager.start_visu, '%Y-%m-%dT%H:%M:%S') and x < datetime.strptime(self.satManager.stop_visu, '%Y-%m-%dT%H:%M:%S') for x in (start, end)):
            errorMsg = wx.MessageDialog(self, "Error: Invalid datetime - Start/End must fall within the maximum analysis interval for this SIFT visualization.", "Error", style=wx.OK|wx.CENTRE)
            errorMsg.ShowModal()
            return

        #Make sure all values are numerical
        if not all((self.check_numeric(llLat), self.check_numeric(llLon), self.check_numeric(urLat), self.check_numeric(urLon))):
            self.create_error_dialog()
            return
        

        llLat = float(llLat)
        llLon = float(llLon)
        urLat = float(urLat)
        urLon = float(urLon)

        #Make sure all values fit within the ranges we want
        if not all(x >= -180 and x < 180 for x in (llLon, urLon)) and not all(x > -90 and x < 90 for x in (llLat, urLat)):
            self.create_error_dialog()
            return

        #Call function to find overpass info in satManager
        try:
            self.satManager.generateOverpassReport(llLon, urLon, llLat, urLat, start, end, self.cbPoints.GetValue())
            successMsg = wx.MessageDialog(self, "CYGNSS overpass report generated.", "CYGNSS Overpass Report", style=wx.OK|wx.CENTRE)
            successMsg.ShowModal()
        except Exception, e:
            logging.error(str(e))
            failureMsg = wx.MessageDialog(self, "CYGNSS overpass report failed.", "CYGNSS Overpass Report", style=wx.OK|wx.CENTRE)
            failureMsg.ShowModal()

        self.Close()



    def on_close(self, event):
        self.MakeModal(False)
        self.Destroy()
        #Ensure that parent window receives focus on close
        wx.GetApp().TopWindow.SetFocus()

    def check_numeric(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False


    def create_error_dialog(self):
        errorString = "Error:  All values must be numeric and between [-180, 180) for Longitude and [-90, 90) for Latitude."
        errorMsg = wx.MessageDialog(self, errorString, "Error", style=wx.OK|wx.CENTRE)

        errorMsg.ShowModal()

if __name__ == '__main__':
    app = wx.App(False)
    cdb.updateColourDB() #Adds additional required colors to db
    if sys.platform.startswith('win'):
        myappid = 'UniversityOfMichigan.CYGNSS-SIFT' # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app.frame = CygnssFrame()
    app.frame.Show()
    app.MainLoop()
