"""
Data analysis for plate reader experimental results of propagating bacterial waves

Aleksandra Sobieska, School of Physics at University of Edinburgh, 2018

The data show optical density values of propagating bacteria on a 384-well microtiter plate in time.
The plate had melted walls between wells to create 16 rows for bacterial propagation.
Two plate readers were used with different recording formats.

The script allows to analyse and visualise the data in various ways, seperated into several modes:

Mode 1
Plots the bacterial wavefront for one row from multiple plates (position vs time).

Mode 2
Calculates correlation coefficients for all rows between 2 plates.

Mode 3
Plots velocity for each row for various plates
The accepted format of row velocity files - a .txt file where each line is for only one speed value for a given row.

Mode 4
Determines the speed of bacterial wavefront for selected start and end time stamps for a given row

"""
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os.path, time
import re

#determine time difference in OD measurements taken by the newer plate reader
def TimeInterval(next_file):
    first_file = str(1).zfill(5) + ".csv"
    second_file = str(next_file).zfill(5) + ".csv"
    file_time1 = time.gmtime(os.path.getmtime(first_file))
    file_time2 = time.gmtime(os.path.getmtime(second_file))
    #convert the tuple date to string that first fmt format
    date1 = str(file_time1[0]) + "-" + str(file_time1[1]) + "-" + str(file_time1[2]) + " " \
    + str(file_time1[3]) + ":" + str(file_time1[4]) + ":" + str(file_time1[5])
    date2 = str(file_time2[0]) + "-" + str(file_time2[1]) + "-" + str(file_time2[2]) + " " \
    + str(file_time2[3]) + ":" + str(file_time2[4]) + ":" + str(file_time2[5])

    fmt = '%Y-%m-%d %H:%M:%S'
    date1 = datetime.strptime(date1, fmt)
    date2 = datetime.strptime(date2, fmt)

    if date1 > date2:
        td = date1 - date2
    else:
        td = date2 - date1
    return td.total_seconds()/3600.0 #returns time in hours

def findFirstDBFFile(): #determine the .dbf file with the earliest measurement
    directory = os.listdir(".")
    all_files = []
    for x in directory:
        if x.endswith("dbf"):
            all_files.append(int(x.replace(".dbf", "")))
    return min(all_files)

#determine time difference in OD measurements taken by Fluostar
def FluostarTimeInterval(abs_start, next_file):
    first_file = str(findFirstDBFFile()) + ".dbf" #determine the file with the earliest measurement
    second_file = str(int(findFirstDBFFile()) + next_file - 1) + ".dbf"
    file_time1 = time.gmtime(os.path.getmtime(first_file))
    file_time2 = time.gmtime(os.path.getmtime(second_file))
    #convert the tuple date to string that first fmt format
    date1 = str(file_time1[0]) + "-" + str(file_time1[1]) + "-" + str(file_time1[2]) + " " \
    + str(file_time1[3]) + ":" + str(file_time1[4]) + ":" + str(file_time1[5])
    date2 = str(file_time2[0]) + "-" + str(file_time2[1]) + "-" + str(file_time2[2]) + " " \
    + str(file_time2[3]) + ":" + str(file_time2[4]) + ":" + str(file_time2[5])

    fmt = '%Y-%m-%d %H:%M:%S'
    date1 = datetime.strptime(date1, fmt)
    date2 = datetime.strptime(date2, fmt)

    if date1 > date2:
        td = date1 - date2
    else:
        td = date2 - date1
    return td.total_seconds()/3600.0 #in hours

#determine the position of bacterial wavefront in time from data taken by Fluostar plate reader
def Fluostar(start, end, r):
    #the OD threshold is the assumed optical density where bacterial wavefront exists
    well = []
    time = []
    OD_positions = []
    all_rows = []

    count = int(start)
    while (count <= end): #throughout all timesteps

        filename = str(count) + str(".csv")
        bacteria_ts = open(filename, "r") #filename should be in double quotes (string)

        lines = bacteria_ts.readlines()
        lines = lines[3:] #removes the last line and the header from the .csv file

        for k in range(0, len(lines)):
            well_measurements = lines[k].split(",")
            OD_positions.append(well_measurements[3])

            if well_measurements[0].endswith("24"): #marks the end of the row
                #convert the list into np array and its elements to floats
                OD_positions = np.array(OD_positions, dtype=np.float32)
                OD_positions = np.fliplr([OD_positions])[0]
                all_rows.append(OD_positions)
                OD_positions = []

        for m in range(0, all_rows[r].size): #in one row
            if all_rows[r][m] <= 0.5: #OD threshold passed or not
                well.append(m)
                time.append(FluostarTimeInterval(start, count))
                break
        count += 2 #avoid luminescence measurements from the folder
        all_rows = []

    bacteria_ts.close()
    return well, time

#determine the position of bacterial wavefront in time from data taken by the newer plate reader
def newPR(start, end, r):
    well = []
    time = []

    count = int(start)
    while (count <= end): #throughout all timesteps
        timestep = str(count).zfill(5)
        filename = timestep + str(".csv")
        bacteria_ts = open(filename, "r") #filename should be in double quotes (string)
        lines = bacteria_ts.readlines()
        lines = lines[:-1] #removes the last line

        #chosen row for the experiment
        OD_positions = lines[r].split(",")
        OD_positions = np.array(OD_positions, dtype=np.float32)
        OD_positions = np.fliplr([OD_positions])[0]
        for m in range(0, OD_positions.size): #in one row
            if OD_positions[m] <= 0.5: #OD threshold passed or not
                well.append(m)
                time.append(TimeInterval(count))
                break
        count1 += 2 #avoid luminescence measurements from the folder

    bacteria_ts.close()
    return well, time

def PRmode(PR, start, end, r): #processing data into lists depending on the plate reader
    #r - row number
    #PR - type of plate reader
    #i - used to distinguish plates by colour
    #well - list with well position of the wavefront
    #time - list with time in hours for corresposning well positions of the wavefront
    if PR == 1:
        well, time = newPR(start, end, r)
    else:
        well, time = Fluostar(start, end, r)

    return well, time

def halveList(well):
    halved_well = []
    for l in range(len(well)):
        if l % 2 == 1:
            halved_well.append(well[l])
    well = halved_well
    return well

def determinePR(plate_directory): #determine which plate reader was used for measurements
    if "Fluostar" in plate_directory:
        PR = 0
    else:
        PR = 1
    return PR

def findLastFile(path_directory, plate_directory): #find the file with the latest measurement
    directory = os.listdir(path_directory + plate_directory)
    all_files = []
    for x in directory:
        if x.endswith("csv"):
            all_files.append(int(x.replace(".csv", "")))
    last_file = max(all_files)
    return last_file

def findFirstCSVFile(path_directory, plate_directory): #determine the .csv file with the earliest measurement
    directory = os.listdir(path_directory + plate_directory)
    all_files = []
    for x in directory:
        if x.endswith("csv"):
            all_files.append(int(x.replace(".csv", "")))
    return min(all_files)

def findSlope(filestart, fileend, row):
    if PR == 1:
        well, time = newPR(filestart, fileend, row)
    else:
        well, time = Fluostar(filestart, fileend, row)

    #uncertainty method taken from https://stackoverflow.com/questions/27634270/how-to-find-error-on-slope-and-intercept-using-numpy-polyfit
    slope, intercept = np.polyfit(time, well, 1, cov=True)
    print "speed: {} +/- {}".format(slope[0], np.sqrt(intercept[0][0]))

def mode1(path_directory):
    #plotting the bacterial wavefront for one row from multiple plates (position vs time)
    plate_amount = int(input("How many plates do you want to compare? "))
    row = int(input("For which row? ")) - 1 #-1 to fit to Python list numbering

    plate_directory = plate_amount*[0]
    PR = plate_amount*[0]
    start = plate_amount*[0]
    end = plate_amount*[0]
    well = plate_amount*[None]
    time = plate_amount*[None]
    labels = []

    for k in range(0, plate_amount):
        plate_directory[k] = raw_input("What is the title of the plate folder " + str(k + 1) + "? ")
        PR[k] = determinePR(plate_directory[k]) #determine which plate reader was used for measurements

        #determine the files with the earliest and latest measurements
        start[k] = int(findFirstCSVFile(path_directory, plate_directory[k]))
        end[k] = int(findLastFile(path_directory, plate_directory[k]))

        os.chdir(path_directory + plate_directory[k])
        well, time = PRmode(PR[k], start[k], end[k], row)

        if k % 5 == 0:
            colour = "r"
        elif k % 5 == 1:
            colour = "b"
        elif k % 5 == 2:
            colour = "c"
        elif k % 5 == 3:
            colour = "y"
        else:
            colour = "g"

        labels.append(plt.scatter(time, well, color = colour, label = plate_directory[k]))


    plt.title("Position of the bacterial wavefront over time")
    plt.xlabel("Time in hours")
    plt.ylabel("Position")
    plt.legend(labels, plate_directory)
    plt.show()

def mode2(path_directory):
    #calculating correlation coefficients for all rows between 2 plates
    plate_directory1 = raw_input("What is the title of the plate folder 1? ")
    PR1 = determinePR(plate_directory1)
    start1 = int(findFirstCSVFile(path_directory, plate_directory1)) #determine the file with the earliest measurement
    end1 = int(findLastFile(path_directory, plate_directory1))

    plate_directory2 = raw_input("What is the title of the plate folder 2? ")
    PR2 = determinePR(plate_directory2)

    start2 = int(findFirstCSVFile(path_directory, plate_directory2)) #determine the file with the earliest measurement
    end2 = int(findLastFile(path_directory, plate_directory2))

    for r in range(0, 16): #through all rows
        os.chdir(path_directory + plate_directory1)
        well1, time1 = PRmode(PR1, start1, end1, r)
        #putting both plate reader experiments approximately at the same time intervals
        if PR1 == 0:
            well1 = halveList(well1)

        os.chdir(path_directory + plate_directory2)
        well2, time2 = PRmode(PR2, start2, end2, r)
        #putting both plate reader experiments approximately at the same time intervals
        if PR2 == 0:
            well2 = halveList(well2)

        if len(well1) > len(well2):
            well1 = well1[:len(well2)]
        else:
            well2 = well2[:len(well1)]

        corr_coeff = np.corrcoef(well1, well2)
        print("Row " + str(r + 1) + ":")

        if math.isnan(corr_coeff[0][1]):
            print "Something is wrong with the plate results."
        else:
            print("The correlation coefficient is " + str(corr_coeff[0][1]) + ".")

        #optional plotting of the correlation
        """
        plt.scatter(well1, well2)
        plt.title("Correlation between plate 1 and plate 2")
        plt.xlabel("Plate 1")
        plt.ylabel("Plate 2")
        plt.show()
        """

def mode3(path_directory): #plotting velocity for each row for various plates
    velocity_folder = raw_input("What is the folder in which you keep files with velocities for each row of a plate? ")
    os.chdir(path_directory + velocity_folder)

    plate_amount = int(input("How many plates do you want to compare? "))
    filename = plate_amount*[0]
    velocities = []
    new_label = []
    rows = np.arange(1, 17)

    for k in range(0, plate_amount): #loop through all plates
        filename[k] = raw_input("What is the title of the plate folder " + str(k + 1) + "? ")
        plate_velocities = open(str(filename[k] + ".txt"), "r") #open the file with row velocities for a given plate
        lines = plate_velocities.readlines()

        for m in range(0, 16):
            velocities.append(lines[m])

        plate_velocities.close()

        velocites = np.array(velocities, dtype=np.float32)

        #pick a different colour for each plate
        if k % 7 == 0:
            colour = "r"
        elif k % 7 == 1:
            colour = "b"
        elif k % 7 == 2:
            colour = "c"
        elif k % 7 == 3:
            colour = "y"
        elif k % 7 == 4:
            colour = "g"
        elif k % 7 == 5:
            colour = "k"
        else:
            colour = "m"

        new_label.append(plt.scatter(rows, velocities, color = colour, label = filename[k]))

        velocities = []

    plt.title("Velocity in every row")
    plt.xlabel("Row number")
    plt.ylabel("Velocity [well position/hour]")
    plt.legend(new_label, filename)
    plt.show()

    os.chdir(path_directory) #change directory to the original main folder




def mode4(path_directory): #determining speed of bacterial wave for selected start and end time stamps for a given row
    # Change directory
    filename = raw_input("What is the title of the plate folder? ")
    row = raw_input("For which row do you want to determine the velocity? ")
    os.chdir(path_directory + filename)

    start = str(input("Start: "))
    end = int(input("End: "))

    #the speed is determined by finding the slope of position of bacterial wavefront in time
    findSlope(start, end, row)

def main():
    path_directory = raw_input("Path directory for the main folder with all the results:")
    print("Mode options:")
    print("1 - plotting the bacterial wavefront for one row from multiple plates (position vs time)")
    print("2 - calculating correlation coefficients for all rows between 2 plates")
    print("3 - plotting row velocities for all 16 rows for various plates")
    print("4 - determining speed of bacterial wavefront for selected start and end time stamps for a given row")
    mode = int(input("What mode do you want to work in?: "))

    if mode == 1:
        mode1(path_directory)

    elif mode == 2:
        mode2(path_directory)

    elif mode == 3:
        mode3(path_directory)

    elif mode == 4:
        mode4(path_directory)

    else:
        print("There is no such mode.")

if __name__ == "__main__":
    main()
