"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 1
"""
from sys import argv, exit
# Assigning Constant Variables
CIRCUIT='.circuit'
END='.end'
# Validating The Number Of Arguments
if len(argv)!=2:
    print('\nUsage: %s <inputfile>' %argv[0])
    exit()
# Validating The File Name
try:
    # Opening And Reading The File
    with open(argv[1]) as f:
        lines=f.readlines()
        start=-1
        end=-2
        # Locating The Beginning And End Of The Circuit By Checking For .circuit And .end
        for line in lines:
            if CIRCUIT==line[:len(CIRCUIT)]:
                start=lines.index(line)
            elif END==line[:len(END)]:
                end=lines.index(line)
                break
        # Validating The Content In The Netlist i.e, Checking If .circuit And .end Are Placed Correctly
        if start>=end or start<0 or end<0:
            print('Invalid circuit definition')
            exit(0)
        # Traverse The Circuit Definition From Last Element To First Element And Print Each Line With Words In Reverse Order
        while end-1>start:
            '''
            Removing Blank Spaces At The Beginning
            Removing Comments After '#'
            Splitting The String Into A List With Space As Separator
            '''
            line1=lines[end-1].split('#')[0].split()
            # Reversing The Order Of The Contents In The Given List
            line2=reversed(line1)
            # Joining The Contents Of The List Using spaces
            line3=' '.join(line2)
            # Printing The Final Line
            print(line3)
            end-=1
        # Closing The File
        f.close()
# Printing Error Message For A Wrong Filename
except IOError:
    print('Invalid file')
    exit()