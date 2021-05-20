"""
Course: EE2703-Applied Programming Lab
Name: Nihal Gajjala
Roll Number: EE19B044
Assignment 2
"""

from sys import argv, exit
import numpy as np

# Assigning Constant Variables
CIRCUIT='.circuit'
END='.end'
AC='.ac'

# Class For One Port elements
class onePortElements():
    # Constructor
    def __init__(self, line):
        self.line = line
        self.tokens = self.line.split()
        self.name = elementType(self.tokens[0])
        self.fromNode = self.tokens[1]
        self.toNode = self.tokens[2]
        # DC Voltage Source
        if len(self.tokens) == 5:
            self.type = 'dc'
            self.value = float(self.tokens[4])
        # AC Voltage Source
        elif len(self.tokens) == 6:
            self.type = 'ac'
            V = float(self.tokens[4])/2
            phase = float(self.tokens[5])
            real = V*np.cos(phase)
            imaginary = V*np.sin(phase)
            self.value = complex(real, imaginary)
        # RLC Circuit
        else:
            self.type = 'RLC'
            self.value = float(self.tokens[3])

# Function To Get Element Name
def elementType(token):
    if token[0] == 'R':
        return 'Resistor'
    elif token[0] == 'L':
        return 'Inductor'
    elif token[0] == 'C':
        return 'Capacitor'
    elif token[0] == 'V':
        return 'Independent Voltage Source'
    elif token[0] == 'I':
        return 'Independent Current Source'

# Function To Return A Dictionary Of Nodes
# The GND node is assigned value 0
def nodeDictionary(circuitDefination):
    dictionary = {}
    nodes = [onePortElements(line).fromNode for line in circuitDefination]
    nodes.extend([onePortElements(line).toNode for line in circuitDefination])
    i = 1
    nodes = list(set(nodes))
    for node in nodes:
        if node == 'GND':
            dictionary[node] = 0
        else:
            dictionary[node] = i
            i += 1
    return dictionary

# Function To Get Key For A Value In The Dictionary
def Key(dictionary, value):
    for key in dictionary.keys():
        if dictionary[key] == value:
            return key

# Function To Make A Dictionary For Each Component Of A Particular Type Of Element
def makeDictionary(circuitDefination, element):
    elementDictionary = {}
    elementNames = [onePortElements(line).tokens[0] for line in circuitDefination if onePortElements(line).tokens[0][0] == element]
    for i,name in enumerate(elementNames):
        elementDictionary[name] = i
    return elementDictionary

# Function To Get Frequency Of AC Source
def Freq(lines):
    frequency = 0
    for line in lines:
        if line[:len(AC)] == AC:
            frequency = float(line.split()[2])
    return frequency

# Function To Get Number of Nodes And Number Of Voltage Sources
def Nodes_And_Voltage_Sources(circuitDefination):
    voltageIndependent = [i for i in range(len(circuitDefination)) if circuitDefination[i].split()[0][0] == 'V']
    VS = len(voltageIndependent)
    N = len(nodeDictionary(circuitDefination))
    return N, VS

# Function To Find The From Or To Of A Given Node
def node(circuitDefination, nodeKey, dictionary):
    x = [(i,j) for i in range(len(circuitDefination)) for j in range(len(circuitDefination[i].split())) if circuitDefination[i].split()[j] in dictionary.keys() if dictionary[circuitDefination[i].split()[j]] == nodeKey]
    return x

# Function For Matrix M and b For A Given Node
def mainMatrix(circuitDefination, frequency, nodeKey, dictionary, voltageDictionary, inductorDictionary, M, b):
    x = node(circuitDefination,nodeKey,dictionary)
    N,VS = Nodes_And_Voltage_Sources(circuitDefination)
    for i in x:
        element = onePortElements(circuitDefination[i[0]])
        elementName = circuitDefination[i[0]].split()[0]
        # Resistor
        if elementName[0] == 'R':
            if i[1] == 1:
                adjKey = dictionary[element.toNode]
                M[nodeKey,nodeKey] += 1/(element.value)
                M[nodeKey,adjKey] -= 1/(element.value)
            if i[1] == 2 :
                adjKey = dictionary[element.fromNode]
                M[nodeKey,nodeKey] += 1/(element.value)
                M[nodeKey,adjKey] -= 1/(element.value)
        # Capacitor
        if elementName[0] == 'C':
            if i[1] == 1:
                adjKey = dictionary[element.toNode]
                M[nodeKey,nodeKey] += complex(0, 2*np.pi*frequency*(element.value))
                M[nodeKey,adjKey] -= complex(0, 2*np.pi*frequency*(element.value))
            if i[1] == 2 :
                adjKey = dictionary[element.fromNode]
                M[nodeKey,nodeKey] += complex(0, 2*np.pi*frequency*(element.value))
                M[nodeKey,adjKey] -= complex(0, 2*np.pi*frequency*(element.value))
        # Inductor
        if elementName[0] == 'L':
            try:
                if i[1] == 1:
                    adjKey = dictionary[element.toNode]
                    M[nodeKey,nodeKey] -= complex(0,1/(2*np.pi*frequency*element.value))
                    M[nodeKey,adjKey] += complex(0,1/(2*np.pi*frequency*element.value))
                if i[1] == 2 :
                    adjKey = dictionary[element.fromNode]
                    M[nodeKey,nodeKey] -= complex(0,1/(2*np.pi*frequency*element.value))
                    M[nodeKey,adjKey] += complex(0,1/(2*np.pi*frequency*element.value))
            except ZeroDivisionError:
                index = inductorDictionary[elementName]
                if i[1] == 1:
                    M[nodeKey,N+VS+index] += 1 
                    M[N+VS+index,nodeKey] -= 1
                    b[N+VS+index] = 0
                if i[1] == 2:
                    M[nodeKey,N+VS+index] -= 1
                    M[N+VS+index,nodeKey] += 1
                    b[N+VS+index] = 0
        # Independent Voltage Source
        if elementName[0] == 'V':
            index = voltageDictionary[elementName]
            if i[1]== 1:
                M[nodeKey,N+index] += 1
                M[N+index,nodeKey] -= 1
                b[N+index] = element.value
            if i[1] == 2 :
                M[nodeKey,N+index] -= 1
                M[N+index,nodeKey] +=1
                b[N+index] = element.value
        # Independent Current Source
        if elementName[0] == 'I':
            if i[1]== 1:
                b[nodeKey] -= element.value
            if i[1] == 2 :
                b[nodeKey] += element.value

# Function To Remove Comments
def removeComments(circuitDefination):
    newDefination = []
    for line in circuitDefination:
        index = len(line)-1
        for i in range(len(line)):
            if line[i] == '#' or line[i] =='#\n':
                index = i
                break
        newDefination.append(line[:index])
    return newDefination

# Main Function
try:
    # Validating The Number Of Arguments
    if len(argv) != 2:
        print('Invalid number of arguments')
        exit()
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
            exit()

    '''
    Removing Comments
    Creating Dictionaries
    Determining Number Of Nodes And Number Of Voltage Sources
    '''
    circuitDefination = lines[start+1:end]
    circuitDefination = removeComments(circuitDefination)
    dictionary = nodeDictionary(circuitDefination)
    voltageDictionary = makeDictionary(circuitDefination, 'V')
    inductorDictionary = makeDictionary(circuitDefination, 'L')
    N,VS = Nodes_And_Voltage_Sources(circuitDefination)
    frequency = Freq(lines)

    # Initializing Matrix M and b
    dimentions = N+VS
    M = np.zeros((dimentions, dimentions), dtype = complex)
    b = np.zeros(dimentions, dtype = complex)

    '''
    DC == True :- DC Voltage Sources
    DC == False :- AC Voltage Sources
    '''
    DC = False 
    if frequency == 0:
        DC = True
        M = np.zeros((dimentions+len(inductorDictionary), dimentions+len(inductorDictionary)), dtype = complex)
        b = np.zeros(dimentions+len(inductorDictionary), dtype = complex)

    # Updating The Matrix M And b For Each Node
    for i in range(len(dictionary)):
        mainMatrix(circuitDefination, frequency, i, dictionary, voltageDictionary, inductorDictionary, M, b)
    M[0] = 0
    M[0, 0] = 1

    # Printing Matrix M And b
    print('The M matrix is :\n', M)
    print('The b matrix is :\n', b)
    # Printing The Node Dictionary
    print('The node dictionary is :', dictionary)

    try:
        # Solving The Equation Mx = b
        x = np.linalg.solve(M,b)    
    except Exception:
        # If The Incidence Matrix Is Singular
        print('The incidence matrix cannot be inverted as it is singular. Please provide a valid circuit definition')
        exit()

    # Printing Voltage At Each Node
    for i in range(N):
        print("The voltage at node {} is {}".format(Key(dictionary, i), x[i]))

    # Printing Current Through Each Independent Voltage Source
    for j in range(VS):
        print('The current through source {} is {}'.format(Key(voltageDictionary, j), x[N+j]))

except Exception:
    print('Invalid file')