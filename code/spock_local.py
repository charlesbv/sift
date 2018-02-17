# -*- coding: utf-8 -*-
#import paramiko # used to be called in SIFT version that run SpOCK on cygnss-sftp-1.engin.umich.edu. Not necessary if SpOCK is run on the local machine
import getpass
import sys
import os
import logging
#from spock_cygnss_spec_download_from_cygnss_server_to_kiosk_soc import * # used to be called in SIFT version that run SpOCK on cygnss-sftp-1.engin.umich.edu. Not necessary if SpOCK is run on the local machine 
from sys import platform
from spock_cygnss_spec_parallel_for_sift import *
from report_coverage_ground_station_for_sift_parallel import *
# mpi = 'C:\\cygwin64\\bin\\orterun'
# print mpi
# os.system("rm -Rf try")
# os.system(mpi + " -np 1 spock.exe try.txt")
# raise Exception

def spock_local(start_date, end_date):
    # start_date = '2017-08-02T00:00:00'
    # end_date ='2017-08-03T00:00:00'

    path_spock_folder = './'

    log_spock_local = path_spock_folder  + 'log_spock_local_start_' + start_date.replace(":","_") + '_end_' + end_date.replace(":","_") + '.txt'
    if ( os.path.isfile("my_ground_stations.txt") == False ):
        os.system("cp ../input_sift/my_ground_stations.txt ./" + " >> " + log_spock_local)
    if ( os.path.isfile("cygnss_geometry_2016.txt") == False ):
        os.system("cp ../input_sift/cygnss_geometry_2016.txt ./" + " >> " + log_spock_local)


    main_input_file_name = path_spock_folder  + 'spock_spec_start_' + start_date + '_end_' + end_date + '.txt'

    #TODO: Date validity checking
    logging.info('*** Running SpOCK propagation ***\nNote: SpOCK produces all satellite and specular point input for SIFT.  This process may take several minutes')
    # !!! TODO:
    # - make sure cygnss_geometry_2016.txt and my_ground_stations.txt are in path_spock_folder (can do that at compilation of SpOCK like wtih egm coeff txt file, antenna.bin etc) -> IDEA: add a cyg: in Makefile (like egm:) where you cp these two files from src to a directory. this directoy could always be the same, for example in input_sift subfolder of sift (since you'd use the optioncyg: only within sift anyway) -> DONE: both files are in input_sift dubfolder of sift. this requires to make Makefile with optioncyg (the same way you call option egm and spec). But call option cyg in make of makeall.sh ONLY if running SpOCK within SIFT (so  not if you're using SpOCK independntly of SIFT)
    # - add path to spock.exe in PATH (do that at compilation, like wehn you add the opath to gsl library (you do that by adding a line in .bach_profile)) -> DONE
    # - add path to spock_cygnss_spec_parallel.py in PATH (by simply putting spock_cygnss_spec_parallel.py in same directory as spock.exe) -> NOT NECESSARY because all py scripts are called within Pythonm not as executables
    # - move spock output files necessary to sift the way you do in spock_cygnss_spec_download_from_cygnss_server_to_kiosk_soc (do that after running spock_cygnss_spec_parallel.py)

    os.chdir(path_spock_folder)
    spock_cygnss_spec_parallel_for_sift( start_date , end_date , 'spec')
    #os.system('python spock_cygnss_spec_parallel_for_sift.py ' + start_date + ' ' + end_date +' spec' + " >> " + log_spock_local)
    report_coverage_ground_station_for_sift_parallel(main_input_file_name.split('/')[-1].replace(":", "_"))
    #os.system('python report_coverage_ground_station_for_sift_parallel.py ' + ' ' + main_input_file_name.split('/')[-1].replace(":", "_") + " >> " + log_spock_local)

    # move SpOCK output files to input_sift
    nb_sc = 8
    start_date_no_time = start_date[0:10].replace(":","_")
    filename_gps_tle = path_spock_folder  + 'gps_' + start_date_no_time + '.txt'
    filename_cygnss_tle = path_spock_folder  + 'cygnss_' + start_date_no_time + '.txt'
    geometry_filename =   '../input_sift/cygnss_geometry_2016.txt'
    input_sift_folder = ['sat_positions', 'spec_positions', 'gsCoverage', 'spock_in','storm_forecasts']
    dir_sift_input = os.path.join('..','input_sift',start_date.replace(':','_') + '_end_' + end_date.replace(':','_'))
    if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
        os.system("mkdir " + dir_sift_input   + " >> " + log_spock_local) 
    else:
        os.system("mkdir -p " + dir_sift_input   + " >> " + log_spock_local) 

    if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
        os.system("mkdir log_spock"    + " >> " + log_spock_local) 
    else:
        os.system("mkdir -p log_spock"  + " >> " + log_spock_local) 


    for ifile in range(len(input_sift_folder)): # for mkdir need to have the path as specified by the os (particularly for windows)
        input_sift_folder[ifile] =  os.path.join(dir_sift_input,'input',input_sift_folder[ifile])
    for ifile in range(len(input_sift_folder)):
        if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
            os.system("mkdir " + input_sift_folder[ifile]  + " >> " + log_spock_local ) 
        else:
            os.system("mkdir -p " + input_sift_folder[ifile]  + " >> " + log_spock_local ) 

    input_sift_folder = ['sat_positions', 'spec_positions', 'gsCoverage', 'spock_in', 'spock_in', 'spock_in', 'spock_in','storm_forecasts'] # need to do it again  (see later in this code)
    dir_sift_input = '../input_sift/' + start_date.replace(':','_') + '_end_' + end_date.replace(':','_') # need to do it again (see later in this code)

    for isc in range(nb_sc):
        output_file_path = path_spock_folder  + "spock/spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + "/spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + str(isc+1) + "/"
        output_file_name = "spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_")  + str(isc+1) + ".txt"
        sift_files = [output_file_path + "interpolated_position_LLA_ECEF_" + output_file_name,
                      output_file_path + "specular_" + output_file_name,
                      output_file_path + "coverage/report_all_by_" + output_file_name,
                      main_input_file_name.replace(":","_"),
                      filename_cygnss_tle,
                      filename_gps_tle,
                      geometry_filename,''] 

        nb_file = len(sift_files)-1
        for ifile in range(nb_file):
            if isc == 0:
                input_sift_folder[ifile] =  dir_sift_input + '/input/' + input_sift_folder[ifile] # here need the separater to be '/' for cp (not for mkfir (see before in the code)) (even in windows)
            filename_only = sift_files[ifile].split('/')[-1].replace(":","_")
            if ((ifile < 3) | (isc == 0)):
                #print "XXX\ncp " + sift_files[ifile] + " " + input_sift_folder[ifile]
                os.system("cp " + sift_files[ifile] + " " + input_sift_folder[ifile] + " >> " + log_spock_local)


    os.system("rm -Rf " + filename_cygnss_tle + " " + filename_gps_tle + " " + main_input_file_name.replace(":","_") + " *_DGD.txt *_DSD.txt " + path_spock_folder  + "spock/spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + " >> " + log_spock_local)
    os.system("mv log_" + 'spock_spec_start_' + start_date.replace(":","_") + '_end_' + end_date.replace(":","_") + " log_spock/" + " >> " + log_spock_local)
    os.system("mv log_gps_tle_" + start_date_no_time + ".txt" + " log_spock/" + " >> " + log_spock_local)
    os.system("mv log_cygnss_tle_" + start_date_no_time + ".txt" + " log_spock/" + " >> " + log_spock_local)
    os.system("mv " + log_spock_local + " log_spock/")
    
    return 

