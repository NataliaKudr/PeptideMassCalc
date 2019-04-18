#Dictionaries and variables
#
 
Da_average = {'A': 71.0779, 'R': 156.1857, 'N': 114.1026, 'D': 115.0874, 'C': 103.1429, 'E': 129.114, 'Q': 128.1292, 'G': 57.0513, 'H': 137.1393, 'I': 113.1576, 'L': 113.1576, 'K': 128.1723, 'M': 131.1961, 'F': 147.1739, 'P': 97.1152, 'S': 87.0773, 'T': 101.1039, 'U': 150.0379, 'W': 186.2099, 'Y': 163.1733, 'V': 99.1311, 'X': 110, 'B': 114.595, 'Z': 128.6216, '*': 0}
Da_MI = {'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943, 'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464, 'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963, 'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028, 'T': 101.047679, 'U': 150.95363, 'W': 186.079313, 'Y': 163.06332, 'V': 99.068414, 'X': 110, 'B': 114.534935, 'Z': 128.550586, '*': 0}
Initial = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'U': 0, 'W': 0, 'Y': 0, 'V': 0, 'X': 0, 'B': 0, 'Z': 0, '*': 0};

water_average = 18.0153
water_MI = 18.0106
proton = 1

#import
#

import sys
from operator import itemgetter
import argparse

#arguments
#

parser = argparse.ArgumentParser(description="Welcome to Peptide Mass Analyser! Please, choose any of the optional arguments to suit your query") #conflict_handler='resolve'
parser.add_argument("-in", "--infile", required=True,   
	help="name of file", type=str)
parser.add_argument("-mass", "--masstocharge", default="a",   
	help="mass to charge type, enter mi for Monoisotopic")
parser.add_argument("-mod", "--modification", 
	help="post-translational modification")
parser.add_argument("-n", "--Nterm", action="store_true",   
	help="N terminal")
parser.add_argument("-c", "--Cterm", action="store_true",  
	help="C terminal")
parser.add_argument("-charge", "--charge",   
	help="charge of ion", type=int, default=1)
parser.add_argument("-cle", "--missedcleavage",   
	help="missed cleavages", type=int)
parser.add_argument("-t", "--top", action="store_true",  
	help="largest 10 by mass to charge")
parser.add_argument("-b", "--bottom", action="store_true",
	help="smallest 10 by mass to charge")
parser.add_argument("-out", "--outfile",   
	help="save to a file", type=str)
parser.add_argument("-show", "--show",   
	help="show the file", action="store_true")
args = parser.parse_args()


# main
#
try:
     fileObj = open(args.infile, 'r') #try to open the file
except FileNotFoundError as errorObj: 
    print("File not found", errorObj)

Dict = {} #creating a dictionary that  will de used as a transit dictionary	for calculations of mass to charge 

# depending on the argument input by the user a dictionary containing mono isotopic or average masses is used
if args.masstocharge == "mi": #Monoisotopic mass, 
    Dict = Da_MI;
    water = water_MI
elif args.masstocharge == "a": #average mass
    Dict = Da_average;
    water = water_average
elif args.masstocharge:
   print("Mass to charge type invalid [mi/a]");
   sys.exit(1) #quit
  
if args.charge: #charge input by the user
    	charge = args.charge
else:
    	charge = 1 #default  
    		

#ptm
#accounting for ptms if chosen by the user
#
if args.modification == "c":  #30.00610 + C -> for S-nitrosylation  	   		
	Dict.update({"C" : 133.149})
elif args.modification == "m": #15.0345 + R and K -> Methylation
	Dict.update({"K" : 143.2068, "R" : 171.2202 }) 	
elif args.modification == "p" : #78.97196 + S T Y -> Phosporylation 
	 Dict.update({"S" : 166.04926, "T" : 180.07586, "Y" : 242.14526  }) 
elif args.modification == "o": #165.2107 + M -> Oxidation of methionines  	   		
	Dict.update({"M" : 296.4068})
elif args.modification:
    print("Modification type invalid [[c|m|p|o]]");
    sys.exit(1) #quit
    	           	
            	

#initiate lists for all elements that will be transcribed into the final outcome
#
name = []
peptide = []
masstocharge = []
charge1 = []
cleavage = []
seq = []
values =[]
term =[]
my_list = [] #create a new list(will be used to append and export new data composed of lists above)

line = fileObj.readlines() #read file
string = ''.join(line) #make a string out of a list
sublines = string.splitlines() #split into lines
data = ' '.join(sublines) #make a string out of each line
final = data.split('>') #split on the appearance of ">"
my_line = final[1:] #remove '' before the first line of interest

for i in my_line:
    item = i.split() #split strings into elements
    pep = item[6] #assign element 6 (amino acid sequence) to a variable 
    if pep == '*': #skip a loop if peptide is merely a stop codon "*", this peptide will not be included in the output 
    	continue
	
    cle = item[3] #assign missed cleavage element
    cle1 = ''.join(cle) #split element 3(missed cleavage) in order to
    num_cle = int(cle1[7]) #extract number of cleavages 
    
   
    for amino in pep:
    	if KeyError in Initial.keys(): 
    		print("Invalid Amino Acid letter discovered: J or O"); #print an error if these 2 amino acids are discovered, every other letter is covered by the dictionaries
    		sys.exit(1)
    	else:
    		Initial[amino] += 1 #add 1 on every appearance of an element in the Initial dictionary
          
    newDict = {} #create a new dictionary to store new values for each  amino acid
    for x in Initial:
    	newDict[x] = Initial[x] * Dict[x] #fill newDict with value derived from multiplication of counts of amino acids and their mass value
    	 
    mass = sum(newDict.values()) #sum all the masses of amino acids present in a sequence
    mtc = (mass + water + proton)/ charge #actual mass to charge value required for the outcome data
    mtcdc = round(mtc,4)  #round mass to charge to 4 decimal places

    def app(): #define the function, where each element is appended, this outcome will be used for calling peptides of interest (c_term, n_term, top, bottom, cleavage args)
    	name.append(item[0]);
    	peptide.append(int(item[2]));
    	masstocharge.append(mtcdc);
    	charge1.append(charge);
    	cleavage.append(num_cle);
    	seq.append(item[6])
    	term.append(item[4])
    	return	
    
    final_str = ('{0:<15}\t{1:>3}\t{2:>10}\t{3:<1}\t{4:>1}\t{5:<100}'.format(item[0].strip(), item[2], mtcdc, str(charge), str(num_cle), item[6])) #formating an outcome list
    my_list.append(final_str) #construct a list made out of values derived every time when each line in a file is read
    
    app() #run app() function
    
    newDict = newDict.fromkeys(newDict, 0) #double-check, reset to zero both dictionaries
    Initial = Initial.fromkeys(Initial, 0) #double-check, reset to zero both dictionaries

fileObj.close() #close the file


out_file = '\n'.join(my_list) #split the list into lines
tuple_my_list = zip(name, peptide, masstocharge, charge1, cleavage, seq, term) #this is a zip file containing the same information as my_list but in a more accessible format for searching

if args.show: #show all data, in the format of the output file
	print(out_file)
	
if args.outfile: #creating, writing to a file with the data on demand
	fh = open(args.outfile,'w')
	print(out_file, file = fh)

if args.top: #look for top 10 biggest peptides
	print("largest 10..", sorted(tuple_my_list, key=itemgetter(2), reverse=True)[:10])
	
if args.bottom: #look for top 10 smallest peptides
	print("smallest 10..", sorted(tuple_my_list, key=itemgetter(2))[:10])

if args.Cterm: #show C term peptides only
	#print("C terminus...", [item for item in tuple_my_list if item[5].endswith("*")]) #calling peptide sequences ending on "*" (was an option when no N, I, C abbreviations in the infile were not available)
	print("C terminus..", [item for item in tuple_my_list if item[6] == "C"])
	
if args.Nterm: #show N term peptides only
	#print("N terminus..", [item for item in tuple_my_list if item[1] == 1])	#calling for peptide 1, where peptide 1 stand for the first peptide in sequence, hence N terminus
	#(was an option when no N, I, C abbreviations in the infile were not available)
	print("C terminus..", [item for item in tuple_my_list if item[6] == "N"])
if args.missedcleavage == 1: #for loop for peptides with n missed cleavages
	print("1 missed cleavage...", [item for item in tuple_my_list if item[4] == 1])
elif args.missedcleavage == 2:
	print("2 missed cleavages...", [item for item in tuple_my_list if item[4] == 2])
elif args.missedcleavage == 0:
	print("0 missed cleavages...", [item for item in tuple_my_list if item[4] == 0])
elif args.missedcleavage == 3:
	print("3 missed cleavages...", [item for item in tuple_my_list if item[4] == 3])
elif args.missedcleavage:
	print("0-3 missed clevages allowed")
	
sys.exit() #exit the program