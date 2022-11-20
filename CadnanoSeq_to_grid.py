#!/usr/bin/python3

#THIS PROGRAM DOESNT TAKE INSERTIONS AND SKIPS INTO ACCOUNT
#Because I have never used them and including them would have considerably complicated the algorithm. 

import sys
import json
import numpy
import pandas as pd
import re 
from itertools import islice
from itertools import accumulate

if __name__ == "__main__":
    cadnanofile = sys.argv[1]
    sequencefile = sys.argv[2]
    outputfilename = sys.argv[3]
    
with open(cadnanofile, 'r') as jf: #importing the json file into python
    cadnano = json.load(jf)

strands = cadnano.get('vstrands') #strands is a list of dictionaries with each dictionary corresponding to one helix

#json file is imported as a dictionary and this code gives the number of items
#of which there are two whose keys are 'name' and 'vstrands'

strands = cadnano.get('vstrands') #vstrands is a key containing scaffold & staple for every helix
m = len(strands[0].get('scaf')) #gives number of boxes filled or otherwise for the helices
n = len(strands) #gives the number of helices
grid = [[' ' for i in range(m)] for _ in range(n)] #creating empty grid which I'll fill later on

class helix:
    def __init__(self, number, boxes):
        self.number = number
        self.boxes = boxes
        
    def __repr__(self):
        return f"helix{self.number}"


allhel = [helix(i, strands[i].get('stap')) for i in range(len(strands))]

staplefile = pd.read_csv(sequencefile) #importing staple sequences file into pandas

# **I thought to instantiate each staple as part of a class and then use its attributes to treat its sequence such that I can apporpriately fill the grid.**


class staple:
    def __init__(self, start, end, sequence):
        self.start = start
        self.end = end
        self.sequence = sequence 
        
        self.startinghelix = int(start.split('[')[0])
        self.endinghelix = int(end.split('[')[0])   
        self.nucleotides = [i for i in self.sequence]
        
        result = re.search(r"\[([0-9]+)\]", self.start)#using regex to find the starting box/square number
        self.start_point = int(result.group(1))
        
        result = re.search(r"\[([0-9]+)\]", self.end)#using regex to find the ending box/square number
        self.end_point = int(result.group(1))

    def absolutegap(boxes, start_point, startinghelix):
        #choosing direction with if and else statements
        if [boxes[start_point + 1][0], boxes[start_point + 1][2]]  == [startinghelix, startinghelix]:
            #go forward and find break 
            stp = boxes[start_point:] #creating a list from which I'll start looping to find when the break occurs in staple seq
            gap = [i for i in range(len(stp)) if stp[i][2] != startinghelix][0] + start_point
        else:
            #go backwards
            stp = boxes[:start_point] #this list ends at the starting point because we have to go backwards in the json file
            gap = [i for i in reversed(range(len(stp))) if stp[i][2] != startinghelix][0]
        return gap

    #we are always going from 5' to 3' on a staple which is why I am using 'if stp[i][2] != self.startinghelix'
    #as a check. Because when the break occurs it's gon be indicated with the connecting helix
    #number on the 3' end i.e. the the third element in the array, [5'H, Box, 3'H, Box].     
    
    def crossover(boxes, interval):
        if boxes[interval][2] != -1:
            #it's a crossover, return the connecting helix
            conn_helix = boxes[interval][2]
            return conn_helix
        else:
            #End of oligo; No crossover here
            conn_helix = None
            
            
    def splitting(self):
        self.firstsplit = staple.absolutegap(allhel[self.startinghelix].boxes, self.start_point, self.startinghelix)
        self.conn_helix = staple.crossover(allhel[self.startinghelix].boxes, self.firstsplit)
        self.splits = [self.firstsplit]
        self.conn_helices = []
        
        if self.conn_helix != None:
            self.secondsplit = staple.absolutegap(allhel[self.conn_helix].boxes, self.firstsplit, self.conn_helix)
            self.conn_helix2 = staple.crossover(allhel[self.conn_helix].boxes, self.secondsplit)
            self.splits.append(self.secondsplit)
            self.conn_helices.append(self.conn_helix)
            
            if self.conn_helix2 != None:
                self.thirdsplit = staple.absolutegap(allhel[self.conn_helix2].boxes, self.secondsplit, self.conn_helix2)
                self.conn_helix3 = staple.crossover(allhel[self.conn_helix2].boxes, self.thirdsplit)
                self.splits.append(self.thirdsplit)
                self.conn_helices.append(self.conn_helix2)
                
                if self.conn_helix3 != None:
                    self.fourthsplit = staple.absolutegap(allhel[self.conn_helix3].boxes, self.thirdsplit, self.conn_helix3)
                    self.conn_helix4 = staple.crossover(allhel[self.conn_helix3].boxes, self.fourthsplit)
                    self.splits.append(self.fourthsplit)
                    self.conn_helices.append(self.conn_helix3)
                
                    if self.conn_helix4 != None:
                        self.fifthsplit = staple.absolutegap(allhel[self.conn_helix4].boxes, self.fourthsplit, self.conn_helix4)
                        self.conn_helix5 = staple.crossover(allhel[self.conn_helix4].boxes, self.fifthsplit)
                        self.splits.append(self.fifthsplit)
                        self.conn_helices.append(self.conn_helix4)
                        
                        if self.conn_helix5 != None:
                            self.sixthsplit = staple.absolutegap(allhel[self.conn_helix5].boxes, self.fifthsplit, self.conn_helix5)
                            self.conn_helix6 = staple.crossover(allhel[self.conn_helix5].boxes, self.sixthsplit)
                            self.splits.append(self.sixthsplit)
                            self.conn_helices.append(self.conn_helix5)
                            
                            if self.conn_helix6 != None:
                                self.conn_helices.append(self.conn_helix6)
                
    def __repr__(self):
        return f"Sequence from {self.start} to {self.end}"
    
    def cutandmark(self):
        
        #first, cutting the sequence as per the available splits for each object
        #then find the relevant positions for pasting each cut
       
        
        self.breakpoints = list(abs(numpy.diff(self.splits)) + 1) 
        self.breakpoints.insert(0, abs(self.start_point - self.firstsplit) + 1)
        self.sequencecuts = [self.nucleotides[j - i:j] for i, j in zip(self.breakpoints, accumulate(self.breakpoints))]        
        
        
        self.points = [self.start_point]
        self.points.extend(self.splits)
        self.helices = [self.startinghelix]
        self.helices.extend(self.conn_helices)
        self.markers = list(zip(self.helices, self.points)) #marker is a pair of numbers indicating helix and box for pasting a sequence cut
        
        return self.sequencecuts
    
    def fill_grid(self, grid):  #finally, paste the splitted sequences 
        
        for i in range(len(self.sequencecuts)):
            squares = allhel[self.markers[i][0]].boxes
            if [squares[self.markers[i][1]+1][0], squares[self.markers[i][1]+1][2]]  == [self.markers[i][0], self.markers[i][0]]:
                #paste forward
                for j in range(len(self.sequencecuts[i])):
                    grid[self.markers[i][0]][self.markers[i][1] + j] = self.sequencecuts[i][j]
            else:
                #paste backward
                for j in reversed(range(len(self.sequencecuts[i]))):
                    grid[self.markers[i][0]][self.markers[i][1] - j] = self.sequencecuts[i][j]
                
#GENERAL FUNCTIONS; not associated with any class
                    
def ss_scaffold(scaffold, staple): #to find single stranded regions in origami; dunno how I'll take their sequence though
    a = [0 if x == [-1,-1,-1,-1] else 1 for x in scaffold]
    b = [0 if x == [-1,-1,-1,-1] else 1 for x in staple]
    c = list(zip(a,b))
    d = [i for i in range(len(c)) if c[i][0] != c[i][1]]
    return ",".join(str(x) for x in d)


def mxn(grid):#converting the grid into mxn text shape 
    seq = ''
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            seq += grid[i][j]
        seq += '\n'
    return seq


def mxnfile(seq):#saving the string grid into a text file 
        with open(f'{outputfilename}_mxngridfile.txt', 'w+') as f:
                f.write(seq)


def complement_base(base):
    """Returns the Watson-Crick complement of a base."""
    
    if base in 'Aa':
        return 'T'
    elif base in 'Tt':
        return 'A'
    elif base in 'Gg':
        return 'C'
    elif base in 'Cc':
        return 'G'

def complement(seq):
    """Compute complement of a sequence."""
    
    # Initialize reverse complement
    comp = ''
    
    # Loop through and populate list with reverse complement
    for base in seq:
        if base not in (' ', '\n'):
            comp += complement_base(base)
        elif base == ' ':
            comp += "T"
        elif base == '\n':
            comp += '\n'
        
    return comp

allstap = [staple(i.Start, i.End, i.Sequence) for i in staplefile.itertuples()] #creating instances of
#class staple from rows of the pandas dataframe. 


for i in range(len(allstap)):
    staple.splitting(allstap[i])


for i in range(len(allstap)):
    staple.cutandmark(allstap[i])
    

for i in range(len(allstap)):
    staple.fill_grid(allstap[i], grid)

str_grid = mxn(grid)
comp_str_grid = complement(str_grid)
mxnfile(comp_str_grid)

for i in range(len(strands)):
    nomatch = ss_scaffold(strands[i].get('scaf'),strands[i].get('stap'))
    if len(nomatch) != 0:     
        print(f"Please note that for helix number {i} there are single stranded scaffold sequences in boxes:", nomatch )
    
