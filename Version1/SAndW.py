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
# 06/02/2015 Moved most steps used to calculate PWMs into calculatePWM function in seedAndWobbleModsLim in order to make this code
#	     and the function more legible.

# A NOTE ON DATA STRUCTURES. This version still uses entirely native python objects. Complex data structures are made up of lists of
# dictionaries, dictionaries of lists, and so on. In order to try and provide some clarity in data structures, variable names are 
# written in camelCase, with a trailing _ followed by the overall structure, using D for dictionary, L for list, so something like
# "kmerRanks_DL" is a Dictionary of Lists, while "areaPWM_DDL" is a Dictionary of Dictionary of Lists.


##################################################################
### seed_and_wobble.pl - translated to SAndW.py 
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
      #Should add further explanations of parameters, but I am not entirely sure of all of them just yet.
      print "\nSeed and wobble will run with the following parameters:\n"
      for i in range(len(paramNames)):
	print str(i+1)+". " + paramNames[i]+": ", params[i]

      inp = raw_input('\nPlease enter the number of any parameter you would like to change,\nor press enter to continue: ')
      if inp != "":
	inpIndex = int(inp)-1
	if 0 <= inpIndex <=5:
	  print "Please enter new value for ",paramNames[inpIndex]," (Currently = ",params[inpIndex],")"
	  params[inpIndex] = raw_input('New value = ')
    
    spotLength = int(params[0])
    startPosition = int(params[1])
    EscoreCutOff = float(params[2])
    printSeeds = str(params[3])
    printTopSeeds = str(params[4])
    topN = int(params[5])
    
    #Could add more error checking, but not seems like overkill. Also not sure if the wrappers here are strictly necessary, but just in case.
    
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

	
  # n = 0  # Lines 145, 148 and 158 can be added, and line 149 commented out for quicker testing when debugging. 
  line = seedFile.readline()
  
  #while n <=2: ##FOR TESTING!!! # take out when done debug
  while line: # disable while testing for quicker results
    line = line.rstrip()
    seedPatterns_L.append(line)
    
    if (len(line) - line.count(".")) != order:
      print "Number of positions in seed", line, "in", seedPatternFile, "does not agree with order = ", order, "\n"
      seedFile.close()
      sys.exit(0)

    #n = n+1  # take out when done debug
    line = seedFile.readline()
  
  seedFile.close()


#################################################################
### Read in intensities and create array of sequences
#################################################################

  dataMatrix_LL = []	# dataMatrix is a list of lists. Each entry contains a list [sequence intensity, sequence] in descending order of intensity
  
  while True: #read input file with extra error handling
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
						   #must use float wrapper, else string comparison is used in line 189
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
  topKmerAreas_D = {}
  topKmerMedians_D = {}
  topKmerZScores_D = {}
  keepFract = 0.5


  for tally in range(len(seedPatterns_L)):
    
    if printSeeds == "yes" or tally == 0:
    
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

    kmerRanks_DL = {}   
    kmerIntensities_DL = {}
    #kmerAreas_D = {}   #Not needed - later call to getTrunctuatedAreas
                        # creates and populates this dictionary


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

        
        # the following code reads all patterns which match the spaced seed, from startPos to end of spaced seed, 
        # from the dataMatrix. Every seen pattern is added to observedPatterns_D
        # If De Bruijns algorithm is properly applied, why not just add all possible patterns which
        # match the seed, rather than bother reading through the file???
        
        currentString = dataMatrix_LL[spotNum][1][i-1:i+spacedSeedWidth-1]
                
        if spacedSeed != revSpacedSeed:   #Check for palindrome
	  fwdElement = currentString[0]		#by setup, the first element of a seed will not be 0
	  revElement = currentString[0]
	  
	  for i in range(1, spacedSeedWidth):
	    fwdElement += (currentString[i] if spacedSeed[i] == '1' else '.')
	    revElement += (currentString[i] if revSpacedSeed[i] == '1' else '.')
	  
	  rcFwdElement = rc(fwdElement)
	  rcRevEl = rc(revElement)
        
	  observedPatterns_D[min(fwdElement,rcFwdElement)] = 0
	  observedPatterns_D[min(revElement,rcRevEl)] = 0
        
        else:	# if it is palindrome, fwd and reverse are the same, don't need to calculate both.
	  fwdElement = currentString[0]
	  
	  for i in range(1, spacedSeedWidth):
	    fwdElement += (currentString[i] if spacedSeed[i] == '1' else '.')
	  
	  rcFwdElement = rc(fwdElement)
        
	  observedPatterns_D[min(fwdElement,rcFwdElement)] = 0
        # VERSION 1.2. [Earlier versions 1.0 and 1.1 not included in this collection]
        # COMPARE TO VERSION 1.1. Version 1.2 much faster -10%  total time. Improvement is due to the use of ''.join method with list
	 #comprehension, which is slower over the number of iterations reuqired here than using for loops as above.

      for i in observedPatterns_D:
	
	if i not in kmerRanks_DL:
	  kmerRanks_DL[i] = []
	  kmerIntensities_DL[i] = []

	kmerRanks_DL[i].append(spotNum+1)  #Again, the +1 is to match with Perl.
	kmerIntensities_DL[i].append(dataMatrix_LL[spotNum][0]) #but here, the list is 0 indexed, whereas it is 1 indexed in perl, so +1 not necessary 

      observedPatterns_D = {}   #clear for next run
      if spotNum % 1000 == 0: print "Spot Number:", spotNum
    
    kmerAreas_D = getKmerTruncArea(kmerRanks_DL, numSpotsInList, keepFract)
    ## This creates and populates kmerAreas_D dictionary!

  ###################################################################
  ## part 2. Calculate median of all log (median), median absolute deviation
  ###################################################################

    logMedian_L = []
    bitSeed = spacedSeed.replace(".","0")[::-1]
    decBitSeed = int(bitSeed,2)   # Perl uses oct() to convert bin to decimal, here we use int(,2)
    
    ### THE FOLLOWING COMMENTED OUT CODE IS UNNECESSARRY.  [Line 326 - 331]  
    # The code using the gapConv2Lets function creates all kmers which fit the pattern of the current seed.
    # then, it checks in line 333 IF [now changed to FOR] the word is in the kmersRaw dictionary.
    # i.e, we are looking for the intersection between KMers Raw dictionary and all Kmers of the shape of the spaced seed.
    # This intersection is simply all kmers in kmersRaw. I am not sure if there is a reason for doing this in PERL, but it is 
    # unnecessary here. However, there may be reasons why it is necessary later - see comments at line 349
    # create ALL gapped patterns of same width, order and pattern as spaced seed.    
    
    #for k in range(0, 4**order):      
    #  word = gapConv2Lets(k, order+1, decBitSeed, spacedSeedWidth) #slight error in gapConv2Letters. Was easier to fix here for now.
    #  revComp = rc(word)

     # if word <= revComp or spacedSeed != revSpacedSeed:  
      #  if word > revComp: word = revComp    
   
    for word in kmerIntensities_DL:
      medianIntensity = median(kmerIntensities_DL[word])
      myLog = log(medianIntensity)
      logMedian_L.append(myLog)
        
    medianLogMedian = median(logMedian_L)
    
    deviations_L = [abs(k - medianLogMedian) for k in logMedian_L]

    medianAbsDev = median(deviations_L)


  #############################################################
  ## Part 3. Print: Word / Enrichment / Median Intensity / Z-Score
  #############################################################
    
    #Here, we do use the gapConv2Lets function, to generate all kmers of the required shape. [see comments starting at line 318 above]
    # WHY?? 1) Because, for output we need the kmers in Alphabetical order, and 2) because, if for some reason a kmer in the total list 
    # has not come up in the dataMatrix, we can show that its area, etc. are N/A.
    # No. 2) should not occur, if the Data Matrix has been set up properly, and the seeds are all within the established prameters.
    # However, it may be useful to have this for testing? I am not convinced.Furthermore, having a line saying N/A, is not much more 
    # useful than not having the line at all.
    # This leaves reason 1). I haven't tested yet, but I think that simply sorting the keys in kmerIntensities_DL will probably be a lot faster
    # than the current process. In other words, this will probably be changed soon.
    
    for k in range(0, 4**order):  ### Check if sorting is faster than working through all seeds...
      word = gapConv2Lets(k, order+1, decBitSeed, spacedSeedWidth) 
      revComp = rc(word)

      if word <= revComp or spacedSeed != revSpacedSeed:
        if printSeeds == "yes" or tally == 0: outFile.write(word+"\t"+revComp+"\t")

        if word > revComp: word = revComp   
                                           
        if word in kmerIntensities_DL:
          medianIntensity = median(kmerIntensities_DL[word])
          myLog = log(medianIntensity)
          zScore = (myLog - medianLogMedian) / (1.4826 * medianAbsDev)   

        
          if printSeeds == "yes" or tally == 0:
            outFile.write( "%.5f\t%.2f\t%.4f\n" %(kmerAreas_D[word], medianIntensity, zScore))

          if kmerAreas_D[word] > EscoreCutOff: 
            topKmerAreas_D[word] = kmerAreas_D[word]
            topKmerMedians_D[word] = medianIntensity
            topKmerZScores_D[word] = zScore

        else:
          if printSeeds == "yes" or tally == 0: 
            outFile.write("NA\tNA\tNA\n")

    outFile.close()
    
    kmerRanks_DL = {}   
    kmerIntensities_DL = {}
    kmerAreas_D = {}


#################################################
# Find top N seeds (adapted from A. Philippakis)
#################################################


  print "Finding top", topN, "seeds \n"
  topNels_LD = []     #List of Dictionarys for Top N elements (dic contains element and value) 
  elVals_LD = []      #List of dictionarys containing the values for all elements 
                     
  topNcounter = 0

  for key in topKmerAreas_D:
    tempDic = {"element" : key, 
               "value" : topKmerAreas_D[key],
               "median" : topKmerMedians_D[key],
               "zscore" : topKmerZScores_D[key]}
    elVals_LD.append(tempDic)
    topNcounter += 1

  if topNcounter < topN: topN = topNcounter  #if there are less keys then the required topN,
                                             # we set topN to the number of keys that there are.
  elVals_LD  = sorted(elVals_LD, key = itemgetter("value"), reverse = True) # itemgetter from operator module

  for N in range(topN):
      tempDic = {"element" : elVals_LD[N]["element"],"value" : elVals_LD[N]["value"]}
      topNels_LD.append(tempDic)

  if printSeeds == "yes":
    
    outputFile2 = outputPrefix+"_"+str(order)+"mers_top_enrichment.txt"
    try:
      outFile2 = open( outputFile2 , 'w')
    except IOError:
      print 'Cannot create top k-mers output file.\n'
      sys.exit(0)
    
    topLabel = str(order)+"-mer"
    outFile2.write(topLabel+"\t"+topLabel+"\tE-score\tMedian\tZ-score\n")
    
    for N in elVals_LD:
      rcEl = rc(N["element"])
      outFile2.write("%s\t%s\t%.5f\t%.2f\t%.4f\n" %(N["element"],rcEl, #continued
                    N["value"],N["median"],N["zscore"]))

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

  areaPWM_DDL = {}
  scaleFactor = log(10000)

  outputFile3 = outputPrefix+"_"+str(order)+"mers_pwm.txt"
  
  try:
    outFile3 = open( outputFile3 , 'w')
  except IOError:
    print 'Cannot create PWM output file.\n'
    sys.exit(0)
  
  
  for z in range(len(topNels_LD)):
    ranking = z + 1
  
    print "Currently on element ranked:\t", ranking ,"\n"
    seed = topNels_LD[z]["element"]
    topEscore = topNels_LD[z]["value"]
    print seed, "\t", topEscore ,"\n"

    outFile3.write(str(ranking) +"\t" + seed + "\t" + str(topKmerAreas_D[seed]) +"\n\n")
    
    #have encapsulated the calculation of the PWM in a function in seedAndWobbleModsLim, with full explanation, rather than leave it here.
    areaPWM_DDL[seed] = calculatePWM(seed, kmerRanks_DL, intensityFile, spotLength, startPosition, scaleFactor,totalPatterns_L)

    outFile3.write('Enrichment score matrix\n\n')
  
    for i in bases:   #more sensible than sorting, because we know what they are.
      line = i+":"
      for j in areaPWM_DDL[seed][i]:
      	line = line + "\t" +str(j)
      line = line + "\n"
      outFile3.write(line)

  	# print Energy matrix for enoLOGOS (?)
    outFile3.write("\nEnergy matrix for enoLOGOS\n\n")
  
    line = "PO"
    for i in range(len(areaPWM_DDL[seed]["A"])): 
      line = line + "\t" + str(i+1)
    line = line + "\n"
    outFile3.write(line)

    for k in bases:
      line = k+":"
      for j in areaPWM_DDL[seed][k]:
        logScaled = -j*log(10000)
        line = line + "\t" + str(logScaled)
      line = line + "\n"
      outFile3.write(line)

### The following code was commented out in PERL. 
### Included only in case we might need it, which is unlikely. 
### Not translated, but would be simple enough to do.

  ##  print OUTPUT3 "\nReverse complement matrix for enoLOGOS\n\n";
  ##  print OUTPUT3 "PO";
  ##  for (my $counter=1; $counter<=($#{$areapwm{$seed}{A}}+1); $counter++) {
  ##    print OUTPUT3 "\t$counter";
  ##  }
  ##  print OUTPUT3 "\n";
  ##  foreach $key (sort keys %{$areapwm{$seed}}) {
  ##    my $compkey;
  ##    if ($key eq "A") {$compkey="T";}
  ##    if ($key eq "C") {$compkey="G";}
  ##    if ($key eq "G") {$compkey="C";}
  ##    if ($key eq "T") {$compkey="A";}
  ##    print OUTPUT3 "$compkey:";
  ##    for (my $y=$#{$areapwm{$seed}{$key}}; $y>=0; $y--) {
  ##        my $logscaled = $areapwm{$seed}{$key}[$y]*(-log(10000));
  ##        print OUTPUT3 "\t$logscaled";
  ##    }
  ##    print OUTPUT3 "\n";
  ##  }

    #Print probability matrix
    outFile3.write("\nProbability matrix\n\n")

    for i in bases:
      line = i+":"
      N = 0 

      for j in areaPWM_DDL[seed][i]:
        numerator = exp(log(10000)*j)
        denominator = sum([exp(log(10000)*areaPWM_DDL[seed][k][N]) for k in bases])
      	#could have done this more like the Perl code, but I think my way is more elegant
      	# and uses less typing, should be same amount of processing
      	
        if denominator != 0:
      		probability = numerator/float(denominator)
      	else: 
      		print "Divide by zero error while creating probability matrix"
      		sys.exit(0)
        # Two checks added which I don't think are strictly necessary.
      	# Added float just to ensure that proper division is used.
      	# Also don't think that I need to check if denominator is zero. 
        # But it's only a few lines of code, so just in case...

      	line = line +"\t" + str(probability)
      	N +=1

      line = line + "\n"
      outFile3.write(line)

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
    Resource testing, which can be True or False - Tests time taken and (not implemented yet)  memory usage. 
	Defaults to False"""
    print ""
    sys.exit(0)

    
  main(sys.argv[1:])   
