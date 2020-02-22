"""
Converting .dbf files to .csv in bulk for plate reader Fluostar

Adapted from https://dbfread.readthedocs.io/en/latest/exporting_data.html
"""
import csv
from dbfread import DBF
import os

def mode1(): #converting DBF files to CSV files in bulk
    start = str(input("Start: "))
    end = int(input("End: "))
    i = 1
    count = int(start)
    while (count <= end):
        fileCSV = str(i) + str(".csv")
        fileDBF = str(count) + str(".dbf")
        myFile = open(fileCSV, 'w')

        table = DBF(fileDBF)
        writer = csv.writer(myFile)

        with myFile:
            writer.writerow(table.field_names)
            for record in table:
                writer.writerow(list(record.values()))
        count += 1
        i += 1

def mode2(): #renumbering DBF files from 1
    i = 1
    for filename in os.listdir("."):
        if filename.endswith("dbf"):
            new_name = str(i) + ".dbf"
            os.rename(filename, new_name)
            i += 1

def main():
    plate_directory = raw_input("What is the title of the plate folder? ")
    path_directory = raw_input("In what path directory do you want to save the plate folder? ")
    os.chdir(path_directory + plate_directory)

    print("Mode options:")
    print("1 - converting DBF files to CSV files in bulk")
    print("2 - renumbering DBF files from 1")
    mode = int(input("What mode do you want to work in?: "))

    if mode == 1:
        mode1()

    elif mode == 2:
        mode2()

if __name__ == "__main__":
    main()
