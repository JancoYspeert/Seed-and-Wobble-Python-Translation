#! usr/bin/python
#SeedAndWobbleModsLim.py
#limited seed and wobble module containing only those functions
#necessary for main seed and wobble program to work

#
# (c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
#

# Python coding by J.Yspeert, 2014
# Original coding as above and noted in the comments before each function.
# abbreviations used in my code:
# trunc: trunctuated
# enrich: enrichment
# Rr: reranked


# change log:
# 16/01/2015. Debug. Sill no producing correct output as per perl script. issues with wobble modules
# 30/01/2015 Debug finally finished. Producing correct output as per perl script. Documentation needs to be updated, various OOP and UI and efficiency tests
# 	     and code understandability improvements can be made.
# 05/02/2015 Added calculatePWM function here - moved it out of SAndW.py [main program] to make both more legible.

### IMPORTS ###
from string import maketrans
import re
import sys
from math import exp, log
from SWclasses import PWM

# useful global variables
rcTrTable = maketrans("acgtACGT","tgcaTGCA")
seedTrTable = maketrans("acgtACGT","11111111")
bases = ["A","C","G","T"]


##########################################################
# 1. rc: Takes a string and returns its reverse complement
##########################################################
def rc(kmer):

	revComp = kmer.translate(rcTrTable)[::-1]
	return revComp


##############################################################################
# 2. getKmerTruncArea: 
# Program that takes all k-mers and their ranks and, for each 
# k-mer, calculates the "area" statistic after dropping a certain fraction from
# the foreground and background.
# Recall that  area = (1/(B+F))(averagebackgroundrank-averageforegroundrank).
##############################################################################
# Perl code written by A. Philippakis, 2006
##############################################################################

def findKmerTruncArea(kmerListsDic, spots, keep): 	# Name + slight funtion Change for OOP. May change further
    

	for key in kmerListsDic:
	  kmerListsDic[key].setArea(truncEnrichU1Array(kmerListsDic[key].ranks, spots, keep))
    

##################################################################################
# 3. truncEnrichU1Array: 
# Function to compute the area statistic in the case of one array, 
# but with the proviso that a fraction of the bottom of the foreground and background 
# will be dropped. This is based on the idea that there will be outliers in the bottom that
# one wants to avoid.  Inputs are a list (not necessarily sorted) of the foreground ranks,
# the number of data points, and the fraction of the foreground/background to consider.
##################################################################################
# Perl code written by A. Philippakis, 2006
##################################################################################


def truncEnrichU1Array(FGList, numPoints, keep):   
  
  # FG = ForeGround, BG = BackGround
  FGList = sorted(FGList)
  origFGSize = len(FGList)
  
  BGSize = int((numPoints-origFGSize) * keep)
  FGSize = int(origFGSize * keep)   		#both after trunctuation 
  
  rankSum = 0
  lastRank = BGSize
  

  for i in range(FGSize):
    
    if (FGList[i]-i) > BGSize:
      rankSum += (lastRank + 1)
      lastRank += 1
      
    else: 
      rankSum += FGList[i]
      lastRank += 1
    
  
  if FGSize > 0 and BGSize != 0:
    result = ((FGSize*(FGSize+1)/2.0) + (FGSize * BGSize/2.0) - rankSum) / float(FGSize*BGSize)
  else:
    if BGSize == 0: print "Division by zero error in truncEnrichU1Array, data may be compromised"
    result = 0
  
  return result


######################################################################
# 4. gapConv2Lets:
# Very similar to convert to letters, but also takes a gapped pattern
# to generate the corresponding word
######################################################################
# convert_to_letters and presumably gapped_convert_to_letters... Perl Code written by A. Philippakis, 2006
#####################################################################################
# Probably not the way I would have done this, but it works. May change it later.

def gapConv2Lets(num, k, seed, width):		# Gapped Convert to Letters
	
	mystr = ""
	sph = 0 	#seed position holder
	domain = [2*(l-1) for l in range(1,k)]
	domain.reverse()

	for i in domain:
		while ((1<<sph)&seed) == 0:
			mystr = mystr + "."
			sph += 1
		mystr +=  bases[((num>>i)&0x3)]
		sph+=1

	return mystr


########################################################
# 5.Median:
# Function returns the median of an array of numbers.
########################################################

def median(myVals):

	length = len(myVals)
	isOdd = length % 2
	tempVals = list(myVals)  #don't necessarilly want to change myVals at this stage (reference would be passed if tempVals = myVals) 
	tempVals.sort()
	
	if isOdd == 1:
		center = (length-1)/2
		median = tempVals[center]
	else:
		center = length/2 #Integer division, towards floor 
		median = (tempVals[center-1] + tempVals[center])/(2.0)
							#but must still subtract 1 as length is 1 higher than last index, zero based indexing 

	return median


################################################################
# 6. wobbleSeedRerank:
# Function to build a pwm based on area medians
################################################################
# Perl code written by A. Philippakis, 2006
################################################################

def wobbleSeedRerank(kmerRanks_DL, center, sPWM, arrayFile, spotLength, startPos, keep):
	
	# collapsedRanks_D = {} 		#rather done each time just before we use collapseRanks.
	# Dont think the undef of PWM_DL is needed here... 
	# I think that this is accomplished by line 392 of S+W.py *** CHECK THIS also RECHECK LINE NUMBER WHEN DONE!
	
	collapsedRanks_D = {}
	
	for i in range(len(center)):
		
		
		if center[i] != ".":
			[Avar, Cvar, Gvar, Tvar] = getVariants(center,i)

			if Avar not in kmerRanks_DL:
				addElement(Avar, kmerRanks_DL, arrayFile, spotLength, startPos)
				
			if Cvar not in kmerRanks_DL:
				addElement(Cvar, kmerRanks_DL, arrayFile, spotLength, startPos)

			if Gvar not in kmerRanks_DL:
				addElement(Gvar, kmerRanks_DL, arrayFile, spotLength, startPos)
			
			if Tvar not in kmerRanks_DL:
				addElement(Tvar, kmerRanks_DL, arrayFile, spotLength, startPos)

			numCollapsedObs=collapseRanks(Avar, Cvar, Gvar, Tvar, kmerRanks_DL, collapsedRanks_D)
			# this also populates collapsedRanks_D
			
			Aval=getRrTruncEnrichU(kmerRanks_DL[Avar],collapsedRanks_D,numCollapsedObs, keep)
			Cval=getRrTruncEnrichU(kmerRanks_DL[Cvar],collapsedRanks_D,numCollapsedObs, keep)
			Gval=getRrTruncEnrichU(kmerRanks_DL[Gvar],collapsedRanks_D,numCollapsedObs, keep)
			Tval=getRrTruncEnrichU(kmerRanks_DL[Tvar],collapsedRanks_D,numCollapsedObs, keep)


			sPWM.setVal("A",i,Aval)
			sPWM.setVal("C",i,Cval)
			sPWM.setVal("G",i,Gval)
			sPWM.setVal("T",i,Tval)

		else:
			sPWM.setVal("A",i,0)
			sPWM.setVal("C",i,0)
			sPWM.setVal("G",i,0)
			sPWM.setVal("T",i,0)

###############################################################################
# 7. getVariants:
# Quick function to get all of the variants at a given position for an
# input word;  NOTE THAT THIS RETURNS VARIANTS THAT ARE LOWER IN LEXICOGRAPHIC
# ORDER WITH RESPECT TO REVERSE COMPLEMENTATION.  NOTE THAT THIS ALSO TRIMS
# OFF .'S THAT BEGIN AND END THE STRING
###############################################################################
# Perl code written by A. Philippakis, 2006
###############################################################################

def getVariants(center,wobblePos):
  
  starts = [wobblePos]
  ends = [wobblePos]
  lengthCenter = len(center)
  revCenter = center[::-1]
  
  for i in bases:
    s = center.find(i)
    e = revCenter.find(i) 
    if s >=0:
      starts.append(s)
    if e >= 0:
      ends.append(lengthCenter-e)
      
  startPos = min(starts)
    
  stopPos = max(ends)

  Avariant = center[startPos:wobblePos]+"A"+center[wobblePos+1:stopPos]
  Cvariant = center[startPos:wobblePos]+"C"+center[wobblePos+1:stopPos]
  Gvariant = center[startPos:wobblePos]+"G"+center[wobblePos+1:stopPos]
  Tvariant = center[startPos:wobblePos]+"T"+center[wobblePos+1:stopPos]
  
  Avariant = min(Avariant,rc(Avariant))
  Cvariant = min(Cvariant,rc(Cvariant))
  Gvariant = min(Gvariant,rc(Gvariant))
  Tvariant = min(Tvariant,rc(Tvariant))
  
  return [Avariant, Cvariant, Gvariant, Tvariant]

#########################################################################################
# 8. addElement
# Function takes a k-mer that may have been ignored in the first pass, scans through the
# array, and adds it as needed.
#########################################################################################
# Perl Code written by A. Philippakis, 2006 and
# modified by M. Berger (4/24/07) to discard positions closest to end of spot
#########################################################################################

def addElement(element, kmerRanks_DL, arrayFile, spotLength, startPos):  ## is kmerRanks DL? Check calls
  # Function adds element to kmerRanks_DL

  rank = 1
  rcElement = rc(element)   

  
  try:
    arrayRead = open(arrayFile, "r")
  except IOError:
    pass #better exception handling?
    print 'cannot open', arrayFile
    sys.exit(0)
 
  #the following has been slightly reworked in order to have less checks within the while loop,
  #we especially don't want to have to check whether the key is already in the Dictionary every time.

  key = min(element,rcElement) 

  if key not in kmerRanks_DL:
    kmerRanks_DL[key] = []
  
  line = arrayRead.readline()

  while line:
    temp = line.rstrip().split("\t")
    
    sequence = temp[-1][(startPos-1):spotLength]  
    
    if re.search(element,sequence, re.I) or re.search(rcElement,sequence, re.I): 
      kmerRanks_DL[key].append(rank)
          
    rank += 1
    line = arrayRead.readline()

  arrayRead.close()


###########################################################################
# 9. collapseRanks:
# Function to take ranks of A variant, C variant G variant and T variant
# and collapse them into a hash that stores their rankings with respect to each other
###########################################################################
# Perl code written by A. Philippakis, 2006
###########################################################################


def collapseRanks(Avar, Cvar, Gvar, Tvar, kmerRanks_DL, collapsedRanks_D):
	
	observedRanks_L = []

	for h in [Avar,Cvar,Gvar,Tvar]:
		for i in kmerRanks_DL[h]:
			tester = 1

			if i in observedRanks_L: tester = 0
			if tester:
				observedRanks_L.append(i)

	observedRanks_L.sort()

	num = 1
	for i in observedRanks_L:
		collapsedRanks_D[i] = num
		num += 1

	return len(observedRanks_L)


#################################################################################
# 10. findMinINfoPos
# Function to take a nascent pwm (i.e., the entries are in areas), and find the
# position of minimum information content. Here it converts areas to a boltzman
# using a scaling factor that is passed to the function.  It then computes entropy at 
# each position using a uniform background, and returns the position that is the least info.
# Note that this also takes as input a seed that tells it which positions to ignore.
#################################################################################
# Perl code written by A. Philippakis, 2006
#################################################################################

def findMinInfoPos(sPWM,seed,scaleFactor):


	seedLength = len(seed)

	tempPWM = PWM(seedLength)
	minEntropy = 5
	
	if seedLength != sPWM.length:
		print "PWM not same size as input seed in findMinInfoPos" 
		sys.exit(0)
	
	#removed check that all rows of PWM are the same length... will not happen in this programming.
	
	for i in range(seedLength):
		

		if seed[i].upper() in bases:
			
			denominator = 0
			entropy = 0

			for j in bases:
				tempPWM.setVal(j,i, exp(sPWM.getVal(j,i)*scaleFactor))

				denominator += tempPWM.getVal(j,i)

			for k in bases:
				tempPWM.setVal(k,i,tempPWM.getVal(k,i)/float(denominator)) # A BIT UGLY
				entropy += tempPWM.getVal(k,i) * log(tempPWM.getVal(k,i),2) # log to base of 2, since  log x / log 2 = log base2 of x.
				## !!! MIGHT NEED TO CHECK NOT TAKING LOG of 0!! probably won't happen because seed[i] must be in bases i.e. not "." i.e. not 0 ?

			entropy += 2

			if entropy < minEntropy:
				minEntropy = entropy
				minEntropyPos = i

	return minEntropyPos


############################################################################
# 11. extSeedAllPatternsRr
# Variation on extend_seed_rerank, but this first discards the minimum
#  information position and examines all other positions corresponding to
#  spaced seeds covered in the array design.
# Requires a separate file of gapped patterns specific to that array.
############################################################################
# Perl Code adapted from extend_seed_rerank by M. Berger (4/24/07)
############################################################################

def extSeedAllPatternsRr(kmerRanks_DL, center, minInfoPos, sPWM, arrayFile,
			 spotLength, startPos, keep, Patterns_L):

 	centerSeed = ""
 	collapsedRanks_D = {}
 	

	newCenter = center[0:minInfoPos]+ "." + center[minInfoPos+1:]
 	centerSeed = newCenter.translate(seedTrTable)

 	for i in range(len(center)):

 		if (center[i] == ".") and (i != minInfoPos):

 			centerSeedT = centerSeed[0:i]+"1"+centerSeed[i+1:]
 			querySeed = trimSeed(centerSeedT)
 			revQSeed = querySeed[::-1]
 			seedPresent = False
 			counter = 0

 			for pat in Patterns_L:			
 				if (querySeed == pat) or (revQSeed == pat):
 					seedPresent = True
 					break


 			if seedPresent:

 				[Avar, Cvar, Gvar, Tvar] = getVariants(newCenter, i)

 				if Avar not in kmerRanks_DL:
 					addElement(Avar, kmerRanks_DL, arrayFile, spotLength, startPos)
 				
 				if Cvar not in kmerRanks_DL:
 					addElement(Cvar, kmerRanks_DL, arrayFile, spotLength, startPos)
 					
 				if Gvar not in kmerRanks_DL:
 					addElement(Gvar, kmerRanks_DL, arrayFile, spotLength, startPos)
 					
 				if Tvar not in kmerRanks_DL:
 					addElement(Tvar, kmerRanks_DL, arrayFile, spotLength, startPos)

 				numCollapsedObs = collapseRanks(Avar, Cvar, Gvar, Tvar, kmerRanks_DL, collapsedRanks_D)

 				Aval=getRrTruncEnrichU(kmerRanks_DL[Avar],collapsedRanks_D,numCollapsedObs, keep)
				Cval=getRrTruncEnrichU(kmerRanks_DL[Cvar],collapsedRanks_D,numCollapsedObs, keep)
				Gval=getRrTruncEnrichU(kmerRanks_DL[Gvar],collapsedRanks_D,numCollapsedObs, keep)
				Tval=getRrTruncEnrichU(kmerRanks_DL[Tvar],collapsedRanks_D,numCollapsedObs, keep)
				

				sPWM.setVal("A",i,Aval)
				sPWM.setVal("C",i,Cval)
				sPWM.setVal("G",i,Gval)
				sPWM.setVal("T",i,Tval)




############################################################################################
# 12. getRrTruncEnrichU
# Function to get the truncated area for a given variant AFTER DOING THE RERANKING PROCEDURE
############################################################################################
# Perl code written by A. Philippakis, 2006
############################################################################################

def getRrTruncEnrichU(varRanks_L,collapsedRanks_D,numObs, keep):

	collapsedMedian_L = []

	for i in varRanks_L:
		collapsedMedian_L.append(collapsedRanks_D[i])

	returnVal = truncEnrichU1Array(collapsedMedian_L, numObs, keep)
	return returnVal


##########################################################
# 13. trimSeed
# Takes a spaced seed element (e.g. ..111.1..111..) and 
# trims uninformative positions from end.
##########################################################
# Perl Code written by M. Berger (04/24/07)
##########################################################

def trimSeed(fullSeed):

	start = 0
	end = 0

	for i in range(len(fullSeed)):
		if fullSeed[i] == ".": start += 1
		else: break

	for i in range(len(fullSeed)):
		if fullSeed[-1-i] == ".": end += 1
		else: break

	end = len(fullSeed) - end

	returnStr = fullSeed[start:end]
	return returnStr

##########################################################
# 14. calculatePWMs
# Takes a seed, our data structures kmerRanks_DL and totalPattterns_L, the intensityFile,
# the paramaters spotLength and startPosition and a scaleFactor, defined in the program as log(10000),
# and uses these to construct a PWM for the seed. Further explanation in text.
# This could have been done fairly easily in main functtion, I have extracted it here to add eplanation.
##########################################################
# Originally written in python by J.Yspeert, 03/02/2015, based on seed_and_wobble.pl code as acknowledged at the start of SAndW.py and this file.
##########################################################
def calculatePWM(seed, kmerRanks_DL, intensityFile, spotLength, startPosition, scaleFactor,totalPatterns_L):
  
  #extend seed by 7 possible bases on each side - this seems somewhat arbitrary, as for extSeedAllPatternsRr, based on the total patterns file, 
  #this could be extended up to 9 possible bases on each side and sometimes be a possible pattern. 
  seed ="."*7 + seed + "."*7
  sLength = len(seed)
  
  #Create a PWM containing 0's in all positions, the length of (extended) seed
  seedPWM = PWM(sLength) # Perl can get away without this, but python - as I have used it - needs it. 
  
  #for each position in the exttended seed which has a base (rather than a '.'), wobbleSeedRerank will get all [A,C,G,T] variants for that
  #position, and give each a value based on the tructuated area for that variant: 
  wobbleSeedRerank(kmerRanks_DL,seed, seedPWM, intensityFile, spotLength, startPosition,1)
  # By reference,  this fill up the areaPWM for all positions of the seed which have a base as there value.

  #Now, for each position in the Seed, we check the information content, and return the position with the least information content.
  minInfoPos = findMinInfoPos(seedPWM,seed,scaleFactor)
  #This does not affect the PWM.
  
  #The minInfoPos is set to a '.' All other positions with a '.' are set in turn to a '1', if the pattern created by this change and the bases remaining
  # after removing the minInfoPos form a pattern contained in the list of Total Patterns, the 1 position is changed to each base, and these variations
  # are ranked in a similar manner to wobbleSeedRerank.
  extSeedAllPatternsRr(kmerRanks_DL,seed, minInfoPos, seedPWM, intensityFile, spotLength, startPosition, 1, totalPatterns_L)
  # By reference, this fills up the positions in the PWM not filled in by wobbleSeedRerank
  
  seedPWM.trim()
  
  return seedPWM
