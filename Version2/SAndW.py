#!/usr/bin/python

#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#

#python coding by J.Yspeert
## Structural and top comments still mostly taken from Perl version... may change over time.
# change log:
# 16/01/2015 Debug. Sill no producing correct output as per perl script. the original seed calculaions work. issues with wobble modules
# 30/01/2015 Debug finished. Producing correct output as per perl script. Documentation needs to be updated, various OOP and UI and efficiency tests
# 	     and code understandability improvements can be made.
# 03/02/2015 Added checkParams, quitOnError and resourceTestingFlags. Still not sure about checking max memory usage. 
#	     Some exception handling and other minor code changes.
# 05/02/2015 Moved most steps used to calculate PWMs into calculatePWM function in seedAndWobbleModsLim in order to make this code
#	     and the function more legible.
# 10/02/2015 Object oriented SAndW completed.

# A NOTE ON DATA STRUCTURES. This version mostly uses classes to represent most complex data structures and objects.
# However, there are still complex data structures made up of lists of dictionaries, dictionaries of lists, and so on. In
# order to try and provide some clarity in data structures, variable names are  written in camelCase, with complex objects followed by
# an underscore _  followed by the overall structure, using D for dictionary, L for list, so something like
# "kmerRanks_DL" is a Dictionary of Lists, while "areaPWM_DDL" is a Dictionary of Dictionary of Lists.
# also, lists or dictionaries generally have "List" or "Dic" as the last part of their name, or "_L", or "_D", respectively.
# Due to various revisions, usage is not consistent betweenan ending "List" and "_L" or "Dic" and "_D"  but this should not impede clarity. 

##################################################################
### SAndW.py - python translation of seed_and_wobble.pl 
###
### Takes as input a two-column list of all "combinatorial" probe
###   intensities and sequences.  MUST BE PRE-ORDERED FROM BRIGHTEST
###   TO DIMMEST.
### Also takes file of gapped patterns (e.g., 11111.111) to
###   consider as candidate seeds.
### Also takes file of all gapped patterns (e.g., 111..11.1.11) that
###   are covered evenly on the universal array design.  This may be
###   the same as or longer than the above file of candidate seeds.
###
### Outputs list of all possible k-mers and their corresponding
###   enrichment score, median signal, and Z score (separate file for each seed).
### Also outputs single integrated file for all seeds above a cutoff.
### Also outputs file of top "N" position weight matrices.
###
### M. Berger 05/02/07 (adapted from A. Philippakis)
###   Modified 04/04/08
### 
### translated into Python by J.Yspeert
##################################################################

## imports

from operator import itemgetter
import sys
from math import log, exp
from seedAndWobbleModsLim import *
import time
from SWclasses import *

# Main function... for command line entry point see: if __name__ == __main__ at the end of this file.

# If you want to use SAndW.main() from another program or the python interpretor (not recommended) , the correct input would be a list:
# [intensity filename, order i.e. number of actual bases or 'occupied positions' in each seed, seed pattern filename, total pattern
# filename, output prefix] the listt could also contain three additional, opttional flags at the end each set to True or False:
# Check Parameters - if True will print out the current parameters of the program and prompt the user tto change or accept them, defaults to False
# quit on Error - if False will allow user to attempt different filenames on an I/O error. Defaults to True
# resource testing - if True will measure the time taken, and - NOT IMPLEMENTED YET - memory usage of the program. Defaults to False.


def main(args):
  intensityFile = str(args[0])
  order = int(args[1])
  seedPatternFile = str(args[2])
  totalPatternFile = str(args[3])
  outputPrefix = str(args[4])
  
  #There probably is a more elegant way to do the following...
  checkParams = False
  quitOnError = True
  resourceTesting = False
  
  if len(args) > 5:
    if str(args[5]).upper()[0] == "T":
      checkParams = True
      
  if len(args) > 6:
    if str(args[6]).upper()[0] == "F":
      quitOnError = False
  
  if len(args) > 7:
    if str(args[7]).upper()[0] == "T":
      resourceTesting = True
      
  # in fact the whole preceeding section might not be that Pythonic, but it does the job.
  

##############################################################################
### PARAMETERS


  if checkParams:
    
    inp = "-"
    paramNames = ["Spot Length", "Start Position", "E-Score CutOff", "Print Seeds ('yes' or 'no')", "Print Top Seeds ('yes' or 'no')", "Number of Top Seeds"]
    params = [36,2,0.25,"yes","yes",3]
    
    while inp != "":
      #Should probably add further explanations of parameters
      print "\nSeed and wobble will run with the following parameters:\n"
      for i in range(len(paramNames)):
	print str(i+1)+". " + paramNames[i]+": ", params[i]

      inp = raw_input('\nPlease enter the number of any parameter you would like to change,\nor press enter to continue: ')
      if inp != "":
	inpIndex = int(inp)-1
	if 0 <= inpIndex <= len(params)-1:
	  print "Please enter new value for ",paramNames[inpIndex]," (Currently = ",params[inpIndex],")"
	  params[inpIndex] = raw_input('New value = ')
    
    spotLength = int(params[0])
    startPosition = int(params[1])
    EscoreCutOff = float(params[2])
    printSeeds = str(params[3])
    printTopSeeds = str(params[4])
    topN = int(params[5])
    
    #Could add more error checking, but it seems like overkill. Also not sure if the wrappers are strictly necessary, but just in case.
    
  else:
    spotLength = 36
    startPosition = 2
    EscoreCutOff = 0.25
    printSeeds = "yes"
    printTopSeeds = "yes"
    topN = 3

### Other things which will be useful later... 
  bases = ["A","C","G","T"]
  
  if resourceTesting:
    startTime = time.time()

###########################################################################
### Read in list of gapped patterns to consider for seeds
###########################################################################
  
  seedPatterns_L = []
  
  while True: #read input file with extra error handling
    try:
      seedFile = open(seedPatternFile, 'r')
      break
    except IOError:
      print 'cannot open seed pattern file:', seedPatternFile
      if quitOnError: sys.exit(0)
      
      inp = raw_input('Enter new seed pattern file name, or change permissions  \nof seed pattern file in a different window and press enter, \nor enter Q to quit: ')
      
      if inp.upper() == "Q": sys.exit(0)
      if inp != "": seedPatternFile = inp  

	
  #n = 0  # Lines 161, 164 and 174 can be added, and line 165 commented out for quicker testing when debugging. 
  line = seedFile.readline()
  
  #while n <=2: ##FOR TESTING!!! # take out when done debug
  while line: # replace with above while testing for quicker results
    line = line.rstrip()
    seedPatterns_L.append(line)
    
    if (len(line) - line.count(".")) != order:
      print "Number of positions in seed", line, "in", seedPatternFile, "does not agree with order = ", order, "\n"
      seedFile.close()
      sys.exit(0)

   # n = n+1  # take out when done debug
    line = seedFile.readline()
  
  seedFile.close()


#################################################################
### Read in intensities and create array of sequences
#################################################################

  dataMatrix_LL = []	# dataMatrix is a list of lists. Each entry contains a list [sequence intensity, sequence] in descending order of intensity
  
  while True: #read intensity input file with extra error handling. THIS COULD BE DONE AS A SEPERATE FUNCTION TO ENHANCE MAIN PROGRAM READABILITY
    try:
      intFile = open(intensityFile, 'r')
      break
    except IOError:
      print 'cannot open seed intensity file:', intensityFile
      if quitOnError: sys.exit(0)
      
      inp = raw_input('Enter different intensity file name, or change permissions  \nof intensity file in a different window and press enter, \nor enter Q to quit: ')
      if inp.upper() == "Q": sys.exit(0)
      if inp != "": intensityFile = inp  
  
  spot = 0
  
  text = intFile.readline()
  
  while text:
    line = text.rstrip().split('\t')
    
    dataMatrix_LL.append([float(line[0]),line[1]]) #line[0] is sequence intensity, line[1] is Sequence
						   #must use float wrapper, else string comparison is used in line 206
    if (spot > 0):
      if (dataMatrix_LL[spot][0] > dataMatrix_LL[spot-1][0]):
	print "Error: Probes are not sorted from brightest to dimmest. \n"
	sys.exit(0)
    
    text = intFile.readline()    
    spot += 1

  intFile.close()


###############################################################################
### Part 1. Calculate median intensity and enrichment score and Z score for each 8-mer
###############################################################################
###The following 3 parts are all one long section.
### This is part 1
### RECHECK WHOLE SECTION

  numSpotsInList = len(dataMatrix_LL)
  topKmersDic = KmerValsDic()
  #topKmerAreas_D = {}; #topKmerMedians_D = {};  topKmerZScores_D = {} # these are in topKmersDic
  keepFract = 0.5  # could possibly be a (user-altterable) parameter?


  for tally in range(len(seedPatterns_L)):
    
    if printSeeds == "yes" or tally == 0:   #move this section to later so that all print sections are together? Clarity?
    
      if printSeeds == "yes":
        outputFile1 = outputPrefix+"_"+str(order)+"mers_"+seedPatterns_L[tally]+".txt"
      
      else: outputFile1 = outputPrefix+"_"+str(order)+"mers.txt"
      
      try:
        outFile = open( outputFile1 , 'w')
      except IOError:
        print 'Cannot create Kmer output file', outputFile1 
        sys.exit(0)  #Not bothering with advanced exception handling in output files for now. It can only be an issue of permissions.
      
      topLabel = str(order) + "-mer"
      outFile.write(topLabel+"\t"+topLabel+"\tE-score\tMedian\tZ-score\n")

    kmersRaw = KmerRawDic()
    #kmerRanks_DL = {}; kmerIntensities_DL = {}; kmerAreas_D = {}   #these three are all in kmersRaw

    observedPatterns_D = {}  # Why is this a dictionary, not a list?
                          #because the same element might show up twice in different areas.
                          #i.e. With a list we would have to check each time whether the element is in the list,
                          # or risk adding it multiple times. With a dictionary, set the element to 0,
                          # which will automatically add it if it is not already in the dic, and do nothing if it is in dic. Faster?
 
    spacedSeed = seedPatterns_L[tally]
    revSpacedSeed = spacedSeed[::-1]
    spacedSeedWidth = len(spacedSeed)
    #split_seed and split_rev_seed not needed in python construction

    print "\nCurrently on spaced seed:", spacedSeed 
    

    for spotNum in range(numSpotsInList):   # Starts at 0, not 1 (as in perl)  due to earlier changes in spot/dataMatrix indexing.
                                            # also, python cannot start a list at index 1 
      for i in range(startPosition,spotLength-spacedSeedWidth+2):  # Apparent change in indexing compared to Perl due to range function

        
        # the following code reads in all patterns, following the spaced seed, from startPos to the spot length, 
        # from the dataMatrix. Every seen pattern is added to observedPatterns_D
        # If De Bruijns algorithm is properly applied, why not just add all possible patterns which
        # match the seed, rather than bother reading through the file???
        
        currentString = dataMatrix_LL[spotNum][1][i-1:i+spacedSeedWidth-1]
                
        if spacedSeed != revSpacedSeed:   #Check for palindrome
	  fwdElement = currentString[0]
	  revElement = currentString[0]
	  
	  for i in range(1, spacedSeedWidth):  #these could be done with list comprehensions and ''.join(), but this is much slower over repeated actions.
	    fwdElement += (currentString[i] if spacedSeed[i] == '1' else '.')
	    revElement += (currentString[i] if revSpacedSeed[i] == '1' else '.')
	  
	  rcFwdElement = rc(fwdElement)
	  rcRevEl = rc(revElement)
        
	  observedPatterns_D[min(fwdElement,rcFwdElement)] = 0
	  observedPatterns_D[min(revElement,rcRevEl)] = 0
        
        else:	# if it is palindrome, fwd and reverse are the same
	  fwdElement = currentString[0]
	  
	  for i in range(1, spacedSeedWidth):
	    fwdElement += (currentString[i] if spacedSeed[i] == '1' else '.')
	  
	  rcFwdElement = rc(fwdElement)
        
	  observedPatterns_D[min(fwdElement,rcFwdElement)] = 0


      for i in observedPatterns_D:
	
	if i not in kmersRaw:
	  kmersRaw[i] = [[spotNum+1],[dataMatrix_LL[spotNum][0]]] # if kmer i isn't in the dictionary, add it, with Rank = spotnum+1, intensity from datamatrix
	else:
	  kmersRaw[i].appendVals([spotNum+1,dataMatrix_LL[spotNum][0]]) # else append these values to the respective lists. 

      observedPatterns_D = {}   #clear for next run
      if spotNum % 1000 == 0: print "Spot Number:", spotNum
    
    kmersRaw.calculateAreas(numSpotsInList, keepFract)
    ## Changed from non OOP version. this should fill in the area for each kmer in kmersRaw. 

  ###################################################################
  ## part 2. Calculate median of all log (median), median absolute deviation
  ###################################################################

    logMedian_L = []
    bitSeed = spacedSeed.replace(".","0")[::-1]
    decBitSeed = int(bitSeed,2)   # Perl uses oct() to convert bin to decimal, here we use int(,2)

    ### THE FOLLOWING COMMENTED OUT CODE IS UNNECESSARRY.    
    # The code using the gapConv2Lets function creates all kmers which fit the pattern of the current seed.
    # then, it checks in line 339 IF [now changed to FOR] the word is in the kmersRaw dictionary.
    # i.e, we are looking for the intersection between KMers Raw dictionary and all Kmers of the shape of the spaced seed.
    # This intersection is simply all kmers in kmersRaw. I am not sure if there is a reason for doing this in PERL, but it is 
    # unnecessary here. However, there may be reasons why it is necessary later - see comments at line 355
    
    ### create ALL gapped patterns of same width, order and pattern as spaced seed.    
    #for k in range(0, 4**order):      # There should be a better way than calling function GConv2Lets in two for loops. Will re-examine. 
    #  word = gapConv2Lets(k, order+1, decBitSeed, spacedSeedWidth) #slight error in gapConv2Letters. Was easier to fix here for now.
    #  revComp = rc(word)

    #  if word <= revComp or spacedSeed != revSpacedSeed:  
    #    if word > revComp: word = revComp    
    
    for word in kmersRaw:  #used to be IF and follow from above, same indentation as last if [if word>revComp] above.
      medianIntensity = median(kmersRaw[word].intensities)
      myLog = log(medianIntensity)
      logMedian_L.append(myLog)
        
    medianLogMedian = median(logMedian_L)
    
    deviations_L = [abs(k - medianLogMedian) for k in logMedian_L]

    medianAbsDev = median(deviations_L)


  #############################################################
  ## Part 3. Print: Word / Enrichment / Median Intensity / Z-Score
  #############################################################
  
    #Here, we do use the gapConv2Lets function, to generate all kmers of the required shape. [see comments starting at line 324 above]
    # WHY?? 1) Because, for output we need the kmers in Alphabetical order, and 2) because, if for some reason a kmer in the total list 
    # has not come up in the dataMatrix, we can show that its area, etc. are N/A.
    # No. 2) should not occur, if the Data Matrix has been set up properly, and the seeds are all within the established prameters.
    # However, it may be useful to have this for testing? I am not convinced.Furthermore, having a line saying N/A, is not much more 
    # useful than not having the line at all.
    # This leaves reason 1). I haven't tested yet, but I think that simply sorting the keys in kmersRaw will probably be a lot faster
    # than the current process. In other words, this will probably be changed soon.
    
    for k in range(0, 4**order):  
      word = gapConv2Lets(k, order+1, decBitSeed, spacedSeedWidth) 
      revComp = rc(word)

      if word <= revComp or spacedSeed != revSpacedSeed:
        if printSeeds == "yes" or tally == 0: outFile.write(word+"\t"+revComp+"\t")

        if word > revComp: word = revComp    ## WHY, what's going on here?
                                           ## ALSO WHY DO WE NEED TWO SEPERATE TESTS?
	if word in kmersRaw:
	  medianIntensity = median(kmersRaw[word].intensities)
          myLog = log(medianIntensity)
          zScore = (myLog - medianLogMedian) / (1.4826 * medianAbsDev)   

        
          if printSeeds == "yes" or tally == 0:
            outFile.write( "%.5f\t%.2f\t%.4f\n" %(kmersRaw[word].area, medianIntensity, zScore))

          if kmersRaw[word].area > EscoreCutOff: 
            topKmersDic[word] = [kmersRaw[word].area, medianIntensity, zScore] 

        else:
          if printSeeds == "yes" or tally == 0: 
            outFile.write("NA\tNA\tNA\n")

    outFile.close()
    
    kmersRaw.clear()


#################################################
# Find top N seeds (adapted from A. Philippakis)
#################################################


  print "Finding top", topN, "seeds \n"
                 
  sortedKeysList = topKmersDic.sortedKeysByArea() # I am not sure if it is faster to create a sorted list of dictionaries - as is done in PERL and
					    # the non-OOP SAndW/py, or justt sort the keys and use the dictionary as is. Since we only need
					    # the sorted list for brief printing functions (for now), I think/hope this is better.

  if len(sortedKeysList) < topN: topN = len(sortedKeysList)  #if there are less keys that pass the e-score cut-off than the required topN,
                                             # we set topN to the number of keys that do pass.

  if printSeeds == "yes":
    
    outputFile2 = outputPrefix+"_"+str(order)+"mers_top_enrichment.txt"
    try:
      outFile2 = open( outputFile2 , 'w')
    except IOError:
      print 'Cannot create top k-mers output file.\n'
      sys.exit(0)
    
    topLabel = str(order)+"-mer"
    outFile2.write(topLabel+"\t"+topLabel+"\tE-score\tMedian\tZ-score\n")
    
    for N in sortedKeysList:
      rcEl = rc(N)
      outFile2.write("%s\t%s\t%.5f\t%.2f\t%.4f\n" %(N,rcEl,topKmersDic[N].area,topKmersDic[N].median,topKmersDic[N].zScore))

    outFile2.close()


#########################################################################################
### Read in list of all gapped patterns covered in universal PBM design (for wobble step)
#########################################################################################

  totalPatterns_L = []
  
  while True: #read input file with extra error handling
    try:
      TPFile = open(totalPatternFile, 'r')
      break
    except IOError:
      print 'cannot open total pattern file:', totalPatternFile
      if quitOnError: sys.exit(0)
      
      inp = raw_input('Enter different total pattern file name, or change permissions  \nof total pattern file in a different window and press enter, \nor enter Q to quit: ')
      if inp.upper() == "Q": sys.exit(0)
      if inp != "": totalPatternFile = inp   
  
  
  text = TPFile.readline()
  
  while text:
    text = text.rstrip()
    totalPatterns_L.append(text)
    numInfoPos = len(text) - text.count(".")

    if numInfoPos != order:
      print "Number of positions in seed %s in %s does not agree with order = %s\n" %(text, totalPatternFile, order)
      sys.exit(0)

    text = TPFile.readline()
    
###########################################################################
### Seed-and-Wobble PWM construction
###########################################################################

  scaleFactor = log(10000)

  outputFile3 = outputPrefix+"_"+str(order)+"mers_pwm.txt"
  
  try:
    outFile3 = open( outputFile3 , 'w')
  except IOError:
    print 'Cannot create PWM output file.\n'
    sys.exit(0)
  
  kmerRanks_DL = {} #kmerRanks is not defined at all in the OOP version. By tthis point in the other versions, it is blank anyway
		    #It is useful here, as as we go through the PWMs, we will probably run into the same seeds a few times.
		    #This may change... see notes about very similar seeds.
     
  for z in range(topN):
    ranking = z + 1
  
    print "Currently on element ranked:\t", ranking ,"\n"
    seed = sortedKeysList[z]
    topEscore = topKmersDic[seed].area
    print seed, "\t", topEscore ,"\n"

    outFile3.write(str(ranking) +"\t" + seed + "\t" + str(topKmersDic[seed].area) +"\n\n")
    

    #have encapsulated the calculation of the PWM in a function in seedAndWobbleModsLim, with full explanation, rather than leave it here.
    esmPWM = calculatePWM(seed, kmerRanks_DL, intensityFile, spotLength, startPosition, scaleFactor,totalPatterns_L)

    outFile3.write('Enrichment score matrix\n\n')
  
    esmPWM.writeFun(outFile3)

  	# print Energy matrix for enoLOGOS (?)
    outFile3.write("\nEnergy matrix for enoLOGOS\n\n")
  
    line = "PO"
    for i in range(esmPWM.length): 
      line = line + "\t" + str(i+1)
    line = line + "\n"
    outFile3.write(line)
    
    esmPWM.writeFun(outFile3,-log(10000))

  ## Removed "Reverse complement matrix for enoLOGOS" code, which was commented out in PERL code, and is in otther SAndW.py versions.

    #Print probability matrix
    outFile3.write("\nProbability matrix\n\n")
    
    probMatrixPWM = esmPWM.makeProbPWM()
    probMatrixPWM.writeFun(outFile3)
  
  outFile3.write("\n\n\n")
  outFile3.close
    
  if resourceTesting:
    stopTime = time.time()
    totalTime = stopTime - startTime
    print "Time taken = ", totalTime, " seconds"


#########################################################################################
### Command line entry point
#########################################################################################

if __name__ == "__main__":

# Remember to put any extra needed error checking in here  
  if len(sys.argv) < 6:
    print "\nUsage requires inputs:"
    print "\tFile containing PBM data (sorted by intensities)"
    print "\twidth of k-mers to inspect"
    print "\tFile containing list of candidate seed patterns (e.g., 11111.111)"
    print "\tFile containing list of all covered gapped patterns (e.g., 111..11.1.11)"
    print "\toutput file prefix\n"
    print "For example: \"python SAndW.py TF_combinatorial.txt 8 query_patterns.txt all_patterns.txt output_prefix\"\n"
    print "Additional flags are:"
    print """	
    Check Parameters, which can be True or False - prints and allows for changing of parameters 
    during program operation.
	Defaults to False 
    Quit On Errors, which can be True or False - when True will quit on most errors - generally I/O,
    when False allows user to atempt to enter alternative filenames or change permissions on files if
    an error is encountered. 
	Defaults to True for server-side use."
    Resource testing, which can be True or False - Tests time taken and memory usage. 
	Defaults to False"""
    print ""
    sys.exit(0)

    
  main(sys.argv[1:])   
