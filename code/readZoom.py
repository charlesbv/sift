# -*- coding: utf-8 -*-
def readZoom(fileloc):
	lines = [line.rstrip('\n') for line in open(fileloc)]
	if not lines:
		return

	lines.pop(0) #First line is instructions
	lines.pop(0) #Second line is header
	lines.pop(0) #Third line is example

	names = []
	min_lons = []
	max_lons = []
	min_lats = []
	max_lats = []
	step_lons = []
	step_lats = []

	for line in lines:
		line = line.split(' ')
		names.append(line[0])
		min_lons.append(int(line[1]))
		max_lons.append(int(line[2]))
		min_lats.append(int(line[3]))
		max_lats.append(int(line[4]))
		step_lons.append(int(line[5]))
		step_lats.append(int(line[6]))



	return names, min_lons, max_lons, min_lats, max_lats, step_lons, step_lats
		
