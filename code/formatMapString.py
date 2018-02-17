# -*- coding: utf-8 -*-
import calendar

#Converts a 'linux' UTC string (YYYY-MM-DDTHH:MM:SS) to a more readable format
def formatDisplayTime(inString):
	year = inString[0:4]
	month_digits = inString[5:7]
	month_txt = calendar.month_abbr[int(month_digits)]
	day = inString[8:10]
	time = inString[11:]

	outString =  day + " " + month_txt + " " + year + " " + time + " UTC"

	return outString

def formatWarningTime(stormInfo):
	outString = ""
	first = True
	if(len(stormInfo)):
		for info in stormInfo:
			if not first:
				outString += "\n"
			inString = info[1]
			year = inString[0:4]
			month_digits = inString[4:6]
			month_txt = calendar.month_abbr[int(month_digits)]
			day = inString[6:8]
			hour = inString[8:] + ":00"
			outString += info[0] + ": " + day + " " + month_txt + " " + hour
			first = False
	else:
		outString = "No Storms Available"


	return outString


	
