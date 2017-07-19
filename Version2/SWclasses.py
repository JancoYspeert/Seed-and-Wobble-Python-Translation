#/usr/bin/python

# SWclasses.py
# Classes for use in OOP seed and wobble.

bases = ["A","C","G","T"]
from math import log, exp
# from seedAndWobbleModsLim import truncEnrichU1Array

def truncEnrichU1Array(FGList, numPoints, keep):   #### THIS SHOULD NOT BE HERE! IT IS ONLY HERE TO PREVENT AN INFINIE REGRESS OF IMPORTS
						   ### I will figure this out later.
						   # COULD JUST ADD CLASSES TO seedAndWobbleModsLim
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


class KmerRawDic:
      def __init__(self):
	self.rawDic = {}
	
      def addItem(self, index):
	self.rawDic[index] = KmerRawStats()
      
      def __contains__(self, item):
	return item in self.rawDic
	
      def __setitem__(self, index, value):
	self.rawDic[index] = KmerRawStats(value[0],value[1])
	
      def __getitem__(self, index):
	return self.rawDic[index]
      
      def __iter__(self):
	for item in self.rawDic:
	  yield item
      def clear(self):
	self.rawDic = {}
	
      def calculateAreas(self, spots, keep):
	for key in self.rawDic:
	  self.rawDic[key].area = truncEnrichU1Array(self.rawDic[key].ranks, spots, keep)
	
class KmerRawStats:
    def __init__(self, ranks = [], intensities = [], area = 0):
      self.ranks = ranks
      self.intensities = intensities
      self.area = area
      
    def appendRank(self,rank):
      self.ranks.append(rank)    
    
    def appendInt(self,intensity):
      self.intensities.append(intensity)
      
    def appendVals(self,vals):
      self.ranks.append(vals[0])
      self.intensities.append(vals[1])
    
    def setArea(self, area):
      self.area = area
      
class KmerValsDic:
      def __init__(self):
	self.valsDic = {}
	
      def addItem(self, index):
	self.valsDic[index] = KmerStats()
      
      def __contains__(self, item):
	return item in self.valsDic
	
      def __setitem__(self, index, value):
	self.valsDic[index] = KmerStats(value[0],value[1],value[2])
	
      def __getitem__(self, index):
	return self.valsDic[index]
      
      def __iter__(self):
	for item in self.valsDic:
	  yield item
	  
      def sortedKeysByArea(self):
	keyList = sorted(self.valsDic, key = lambda x: self.valsDic[x].area, reverse = True)
	return keyList
	# There is probably also a better way to do this in the non-OOP SAndW's!]
	
      def clear(self):
	self.valsDic = {}
	
class KmerStats:
  # values: area, median, Zscore
  # There seems to be little reason, at the moment to keep KmerStats and KmerValsDic seperate.
  # However, in future, it may sbe useful if we want to add further stats here.
  
  def __init__(self, area = 0, median = 0, zScore = 0):
    self.area = area
    self.median = median
    self.zScore = zScore

  def setArea(self, area):
    self.area = area
 
  def setMedian(self, median):
    self.median = median
  
  def setZScore(self, zScore):
    self.zScore = zScore

class PWM:		#one could also probably import a matrix class from one of the maths/science modules.
  # Keeping this simple for now, but can add a lot of methods which may be useful for analysis later
  def __init__(self, length):
    self.PWMarray = {k: [0]*length for k in bases}
    self.length = length
    #SHOULD PROBABLY ADD SEED HERE, as itt will be important later, if we compare PWMs constructed using S+W methods.
    
  def setVal(self, base, index, value):
    self.PWMarray[base][index] = value
     # there probably is a way to do this with __setattr__, but it's beyond me at the moment. 
     #Or, one could set up the PWM as just a list of length 4*self.length here, but that's a little bit ugly.
  
  def getVal(self, base, index):
    return self.PWMarray[base][index]
  
  def trim(self):			#removes leading and ttrailing elements which are all zero.
    while [self.PWMarray[k][0] for k in bases] == [0,0,0,0]:
      for l in bases: self.PWMarray[l].pop(0)

    while [self.PWMarray[k][-1] for k in bases] == [0,0,0,0]:
      for l in bases: self.PWMarray[l].pop()
    
    self.length = len(self.PWMarray["A"])
    
  #The following two functions could possibly be combined, but they are already quite messy as it is. Will be improved.  
    
  def writeFun(self,outputHandle,factor =1):	#could overload __repr__ or __str__ functions, but it makes more sense this way, since we are using vaious parameters.
				# outputHandle is the reference to the already opened function file.
				# factor is a factor which each entry in PWM is to be multiplied by.
      for i in bases:
	line = "".join([str(factor * self.PWMarray[i][n])+"\t" for n in range(self.length)])
	outputHandle.write(i+":\t"+line+"\n")

	      
  def makeProbPWM(self):		# creates a probability (normal) PWM from the Enrichment Score Matrix. This is the type we normaly work with
    l = self.length
    probPWM = PWM(l)
    denoms = [0]*l
    for n in range(l):
      denoms[n] = sum([exp(log(10000)*self.PWMarray[k][n]) for k in bases])
    for k in bases:
      for n in range(l):
	probPWM.setVal(k,n,exp(log(10000)*self.PWMarray[k][n])/denoms[n])
    return probPWM
    