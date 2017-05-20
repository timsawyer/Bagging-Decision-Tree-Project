# This Python script was written using Python 2.7.13
# This is the same decision tree that I used in homework one
# but modified to take input from my bagging script

from __future__ import division
from math import log
import scipy.stats
import sys

# python library for loading and reading in arff files
# available from https://pypi.python.org/pypi/liac-arff
# install with: pip install liac-arff
import arff

if len(sys.argv) > 1:
	trainingDataFilename = sys.argv[1]
else:
  print "No training data specified. Exiting."
  exit()

# trainingDataFilename = 'bagging-data/molecular-biology_promoters_train.arff'

class Node:
  def __init__(self, value, name, index, metadata):
    self.value = value # every node is child of classifier attribute - value of the attribute in this branch e.g. Sunny, Overcast, etc.
    self.name = name # name of this attribute e.g. Humidity, Wind
    self.index  = index # index of attribute in data so we know where to find it in each sample
    self.metadata = metadata # metadata for attribute from arff file. 0,1 if boolean
    self.type = "classifier"
    self.children = []

  # define how we convert this object to string so we can print our tree for analysis
  def __str__(self):
    return "{\"type\": \"classifier\", \"index\": \"" + str(self.index) + "\", \"name\":\"" + self.name.encode('utf-8').translate(None, '\\') + "\",\"value\":\"" + self.value.encode('utf-8').translate(None, '\\') + "\",\"children\":[" + reduce(lambda x,y: x + "," + y, map(lambda x: str(x), self.children)) + "]}" 

class Leaf:
  def __init__(self, value, result):
    self.value = value # every node is child of classifier attribute - value of the attribute in this branch e.g. Sunny, Overcast, etc.
    self.result = result # 0, 1 
    self.type = "leaf"
  
  # define how we convert this object to string so we can print our tree for analysis
  def __str__(self):
    return "{\"type\": \"leaf\", \"result\": \"" + str(self.result) + "\", \"value\": \"" + self.value.encode('utf-8').translate(None, '\\') + "\"}"

def getNumNegativeSamples(S):
  """ Returns number of negative samples in a set
  Args:
    S: list of samples
  """
  negSamples = list(filter(lambda x: x[sampleResultIndex] == u'-', S))

  # if no positive samples return 0
  if (len(negSamples) == 0):
    return 0
  else:
    # some samples are fractional, so we must add up sizes of each sample
    numNegSamples = reduce(lambda x, y: x + y , map(lambda x: x[sampleSizeIndex], negSamples))
    return numNegSamples

# return number of positive samples in a set
def getNumPositiveSamples(S):
  """ Returns number of positive samples in a set
  Args:
    S: list of samples
  """
  posSamples = list(filter(lambda x: x[sampleResultIndex] == u'+', S))

  # if no positive samples return 0
  if (len(posSamples) == 0):
    return 0
  else:	
    # some samples are fractional, so we must add up sizes of each sample
    numPosSamples = reduce(lambda x, y: x + y, map(lambda x: x[sampleSizeIndex], posSamples))
    return numPosSamples

def getNumSamples(S):
  """ Returns number of samples in a set
  Args:
    S: list of samples
  """
  if(len(S) == 0):
    return 0
  else:
    # some samples are fractional, so we must add up sizes of each sample
    return reduce(lambda x, y: x + y, map(lambda x: x[sampleSizeIndex], S))

# get the samples from this set that are missing this attribute
def getMissingSamplesOfAttribute(S, attributeIndex):
  return list(filter(lambda x: (x[attributeIndex] == None), S))

# from a set get the samples that have a certain value of an attribute
# e.g. Samples that have an Outlook of Sunny
# for samples that have an unknown(?) value for requested attribute, we calculate fractional amounts to return
def getSubsetWithAttributeValue(S, attributeIndex, attributeValue, numPossibleAttributeValues):
  """ Returns a list
  Args:
    S: list of samples
    attributeIndex: integer that specifies that index of the attribute in the data
    attributeValue: string of the attribute value to search for
    numPossibleAttributeValues: integer that represents total number of possible attribute values
  """
  # get samples that are missing the current best attribute
  S_missing = getMissingSamplesOfAttribute(S, attributeIndex)
  numMissing = len(S_missing)
  numWithAttribute = len(S) - numMissing
  subsetWithAttributeValue = list(filter(lambda x: x[attributeIndex] == attributeValue, S))
  subsetProbability = 0
        
  # add to each subset all the samples that are missing this attribute but
  # assign a fraction to them for each subset depending on the probability that the sample is part of this subset
  # if all samples are missing this attribute, we divide samples evenly among number of possible attribute values
  if (numWithAttribute == 0):
    subsetProbability = 1 / numPossibleAttributeValues
  else:
    subsetProbability = len(subsetWithAttributeValue) / numWithAttribute
  
  # if probability is zero, we don't need to add any
  if (subsetProbability > 0):
    missingAddition = list(map(lambda x: setSizeOfListSamples(x, subsetProbability), S_missing))
    return (subsetWithAttributeValue + missingAddition)
  else:
    return subsetWithAttributeValue

def calcEntropy(S):
  """ Returns a float
  Args:
    S: list of samples
  """
  numSamples = getNumSamples(S)
  if (numSamples == 0):
    return 0

  numNegativeSamples = getNumNegativeSamples(S)
  numPositiveSamples = getNumPositiveSamples(S)
  negativeProportion = numNegativeSamples / numSamples
  positiveProportion = numPositiveSamples / numSamples
  log2neg = log(negativeProportion, 2) if negativeProportion > 0 else 0
  log2pos = log(positiveProportion, 2) if positiveProportion > 0 else 0
  return (-positiveProportion * log2pos) - (negativeProportion * log2neg)

# calculate the information gain of choosing this attribute for current set
def calcGain(S, A, A_index):
  """ Returns a float
  Args:
    S: list of samples
    A: attribute in the form generated from liac-arff library
    A_index: integer that specifies that index of the attribute in the data
  """
  entropyS = calcEntropy(S)
  lengthS = getNumSamples(S)
  entropySumV = 0
  attributeValues = A[1]
  numAttributeValues = len(attributeValues)

  for value in attributeValues:
    S_v = getSubsetWithAttributeValue(S, A_index, value, numAttributeValues)
    lengthS_v = getNumSamples(S_v)
    entropySumV += (lengthS_v / lengthS) * calcEntropy(S_v)
  return entropyS - entropySumV

def passesChiSquareTest(S, attribute, attributeIndex, confidenceLevel):
  """ Returns a boolean
  Args:
    S: list of samples
    A: attribute in the form generated from liac-arff library
    attributeIndex: integer that specifies that index of the attribute in the data
    confidenceLevel: float that represents confidence level needed to reject null hypothesis
  """
  sum = 0
  attributeValues = attribute[1]
  numAttributeValues = len(attributeValues)
  dof = numAttributeValues - 1
  critcalValue = scipy.stats.chi2.isf(1-confidenceLevel, dof)
  numNegativeSamples = getNumNegativeSamples(S)
  numPositiveSamples = getNumPositiveSamples(S)
  positiveSamples = list(filter(lambda x: x[sampleResultIndex] == u'+', S))
  negativeSamples = list(filter(lambda x: x[sampleResultIndex] == u'-', S))

  for value in attributeValues:
    p_samples = getSubsetWithAttributeValue(positiveSamples, attributeIndex, value, numAttributeValues)
    n_samples = getSubsetWithAttributeValue(negativeSamples, attributeIndex, value, numAttributeValues)
    p_i = getNumSamples(p_samples)
    n_i = getNumSamples(n_samples)
    p_prime_i = numPositiveSamples * ((p_i + n_i) / (numPositiveSamples + numNegativeSamples))
    n_prime_i = numNegativeSamples * ((p_i + n_i) / (numPositiveSamples + numNegativeSamples))

    # protection for division by 0. 
    # In this case the numerator of the equation below will be 0 also, so changing to n_prime_i or p_prime_i to 1 doesnt modify outcome
    p_prime_i = 1 if p_prime_i == 0 else p_prime_i
    n_prime_i = 1 if n_prime_i == 0 else n_prime_i

    sum += ((p_i - p_prime_i)**2 / p_prime_i) + ((n_i - n_prime_i)**2 / n_prime_i)

  # if the sum is greater than the critcalValue, then we are confident that we can reject the null hypothesis
  if (sum >= critcalValue):
    return True
  else:
    return False

def chooseBestAttribute(S, availableAttributes, confidenceLevel):
  """ Returns a boolean
  Args:
    S: list of samples
    availableAttributes: list of strings that specify names of attributes that are not in tree yet
    confidenceLevel: float that represents confidence level needed to reject null hypothesis
  """
  if(availableAttributes):
    attributeHG = -1
    attributeHGIndex = -1 # preserve index of attribute, so where know where in data to reference it
    highestGain = -1

    # we always loop through the entire training data attributes
    # this preserves the correct index of the attribute for referencing the data
    # we simply skip calculating gain for an attribute if it has already been used in above tree
    for index, attribute in enumerate(trainingData[u'attributes']):
      if(attribute[0] in availableAttributes):
        # attribute[0] contains name of attribute, attribute[1] contains list of possible attribute values
        # calculate gain for this attribute
        attributeGain = calcGain(S, attribute, index)

        # if higher than any previous gain for another attribute
        # save this attribute as current highest
        if(attributeGain > highestGain):
          highestGain = attributeGain
          attributeHG = attribute
          attributeHGIndex = index

    # after finding attribute with highest gain
    # run chi-square test to see if we want to use it
    if (passesChiSquareTest(S, attributeHG, attributeHGIndex, confidenceLevel)):
      return {
        "metadata": attributeHG,
        "index": attributeHGIndex
      }
    else:
      return
  else:
    return

# recrusive function to build out decision tree
def growTree(S, availableAttributes, parentAttributeValue, confidenceLevel):
  """ Returns a boolean
  Args:
    S: Set to grow tree from
    availableAttributes: attributes that have not been used in parent tree
    parentAttributeValue: value of attribute that this branch is off of. e.g. Sunny, Overcast, Rainy
  """
  if (all(y[sampleResultIndex] == u'-' for y in S)):
    return Leaf(parentAttributeValue, 0)
  
  elif (all(y[sampleResultIndex] == u'+' for y in S)):
    return Leaf(parentAttributeValue, 1)

  else:
    # get best attribute
    bestAttribute = chooseBestAttribute(S, availableAttributes, confidenceLevel)
    if (bestAttribute):
      attributeName = bestAttribute["metadata"][0]
      attributeValues = bestAttribute["metadata"][1]
      numAttributeValues = len(attributeValues)

      # remove this attribute from use in subtree
      subTreeAvailableAttributes = list(filter(lambda x: x != attributeName, availableAttributes))

      # create new node in decision tree
      newNode = Node(parentAttributeValue, attributeName, bestAttribute["index"], bestAttribute["metadata"])

      # split set into subsets based on current best attribute
      for i, value in enumerate(attributeValues):
        subsetWithAttributeValue = getSubsetWithAttributeValue(S, bestAttribute["index"], value, numAttributeValues)

        # for each value of this attribute grow its subtree, if there are any samples for it
        if (len(subsetWithAttributeValue) > 0):
          newNode.children.append(growTree(subsetWithAttributeValue, subTreeAvailableAttributes, value, confidenceLevel))
      return newNode

    else:
      # if we have a mixed group, but no more attributes to choose from (or if no attributes pass chi-squared test)
      # choose true,false depending which has higher percentage in group
      numNegativeSamples = getNumNegativeSamples(S)
      numPositiveSamples = getNumPositiveSamples(S)
      value = 1 if numPositiveSamples > numNegativeSamples else 0
      return Leaf(parentAttributeValue, value)

# create decision tree with confidence level
def createDecisionTree(confidenceLevel):
  # to simply, we only record names of available attributes.
  # All other metadata not needed for just remembering which have already been used
  # x[0] is name of attribute
  # x[1] is list of possible values for that attribute
  # to start all attributes are available
  availableAttributes = list(map(lambda x: x[0], trainingData[u'attributes']))

  # remove first item in list, since it contains the results of sample, i.e. positive or negative
  # therefore is not for use as a classifier
  availableAttributes.pop(0)

  return growTree(sampleSet, availableAttributes, "", confidenceLevel)

# recursive method to traverse decision tree and make prediction for each test sample
def makePrediction(sample, decisionTree):
  if(decisionTree.type == "leaf"):
    return decisionTree.result
  
  else:
    classifierIndex = decisionTree.index
    sampleClassifierValue = sample[classifierIndex]

    for subtree in decisionTree.children:
      if(subtree.value == sampleClassifierValue):
        return makePrediction(sample, subtree)

# function to add default size of sample to each sample
def appendSizeToSample(lst):
  lst.append(1)
  return lst

def setSizeOfListSamples(lst, size):
  lst[sampleSizeIndex] = size
  return lst

# function to use decision tree over test set
def runTestWithDecisionTree(decisionTree, confidenceLevel):
  # now run decision tree against test set and evaluate performance
  numPredictions = len(testData[u'data'])
  numCorrectPredictions = 0

  for sample in testData[u'data']:
    prediction = makePrediction(sample, decisionTree)
    # print str(prediction)
    actualSampleResult = 1 if sample[sampleResultIndex] == u'+' else 0
    if (prediction == actualSampleResult):
      numCorrectPredictions += 1
  print "Percentage of correct predictions using decision tree with confidence level of " + str(100 * confidenceLevel) + "%: " + str(100 * numCorrectPredictions / numPredictions) + "%\n"


# create decision tree, run test set against it, and save it to file
def runTest(confidenceLevel):
  decisionTree = createDecisionTree(confidenceLevel)
  runTestWithDecisionTree(decisionTree, confidenceLevel)
  filename = "DecisionTree_" + str(int(confidenceLevel * 100)) + ".json"
  # saveDecisionTreeToFile(decisionTree, filename)

# load training data and test data
trainingData = arff.load(open(trainingDataFilename, 'rb'))
testData = arff.load(open('bagging-data/molecular-biology_promoters_test.arff', 'rb'))

# index of field that records whether sample is positive or negative
# in this data set, always first value in data
sampleResultIndex = 0

# adding an additional attribute to each sample which is its size
# every sample starts with size 1, but as we run into samples with missing attributes, we create fractional samples
sampleSet = list(map(lambda x: appendSizeToSample(x), trainingData[u'data']))
sampleSizeIndex = len(trainingData[u'attributes']) # index in each sample where to find size of sample

# run test for each confidence level and use it against test set
# this will take a long time if you run all of them
# runTest(.70)

confidenceLevel = .75
decisionTree = createDecisionTree(confidenceLevel)
print str(decisionTree)

# print "Tests complete"




