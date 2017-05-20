# This Python script was written using Python 2.7.13

from __future__ import division
from math import log
import scipy.stats
import random
import copy
import os
import glob
import subprocess
from asyncproc import Process
import time
import json

# python library for loading and reading in arff files
# available from https://pypi.python.org/pypi/liac-arff
# install with: pip install liac-arff
import arff

decisionTreeProcesses = []

# index of field that records whether sample is positive or negative
# in this data set, always first value in data
sampleResultIndex = 0

"""
  Function to remove previously generated replicates
"""
def removeReplicates():
  files = glob.glob('bagging-data/replicates/*')
  for f in files:
    os.remove(f)

"""
  Function to create a bootstrap replicate from the training data
"""
def createBootstrapReplicate(trainingData, filename):
  numTrainingSamples = len(trainingData[u'data'])
  replicateData = copy.deepcopy(trainingData)
  replicateData[u'data'] = []

  for i in range(0, numTrainingSamples):
    randomIndex = random.randint(0, numTrainingSamples-1)
    replicateData[u'data'].append(copy.deepcopy(trainingData[u'data'][randomIndex]))

  outputFile = open(filename, 'wb')
  arff.dump(replicateData, outputFile)
  outputFile.close()


"""
  Function to create a decision tree from a replicate
"""
def createDecisionTree(replicateFilename):
  args = ['python', 'decisionTree.py', replicateFilename]
  decisionTreeProcesses.append(Process(args))

"""
  Function to check if any of our decision tree threads are still running
"""
def areTreesRunning(statusList):
  for treeStatus in statusList:
    if treeStatus:
      return True
  return False

"""
  recursive method to traverse decision tree and make prediction for each test sample
"""
def makePrediction(sample, decisionTree):
  if(decisionTree[u'type'] == u'leaf'):
    result = 1 if decisionTree[u'result'] == u'1' else 0
    return result
  
  else:
    classifierIndex = int(str(decisionTree[u'index']))
    sampleClassifierValue = sample[classifierIndex]

    for subtree in decisionTree[u'children']:
      if(subtree[u'value'] == sampleClassifierValue):
        return makePrediction(sample, subtree)

"""
 Function to use decision trees over test set, and tally votes from each tree
 Majority vote decides prediction
"""
def tallyVotes(votes):
  numPositivites = 0
  numNegatives = 0
  for vote in votes:
    if vote == 0:
      numNegatives += 1
    else:
      numPositivites += 1

  if numPositivites > numNegatives:
    return 1
  else:
    return 0

"""
 Function to use decision trees over test set, and tally votes from each tree
 Majority vote decides prediction
"""
def runTestWithDecisionTrees(decisionTrees):
  # Load test data
  testData = arff.load(open('bagging-data/molecular-biology_promoters_test.arff', 'rb'))

  numPredictions = len(testData[u'data'])
  numCorrectPredictions = 0

  for sample in testData[u'data']:
    votes = []
    for decisionTree in decisionTrees:
      prediction = makePrediction(sample, decisionTree)
      votes.append(prediction)

    actualSampleResult = 1 if sample[sampleResultIndex] == u'+' else 0
    # winningVote = makePrediction(sample, decisionTrees[0])
    winningVote = tallyVotes(votes)
    if (winningVote == actualSampleResult):
      numCorrectPredictions += 1
  print "Percentage of correct predictions using " + str(len(decisionTrees)) + " samplings: " + str(100 * numCorrectPredictions / numPredictions) + "%\n"

# Remove previously generated replicates
removeReplicates()

# Load training data
trainingData = arff.load(open('bagging-data/molecular-biology_promoters_train.arff', 'rb'))

# Set number of samples we're going to do
numSamplings = 20

for i in range(0, numSamplings):
  replicateFilename = 'bagging-data/replicates/replicate_' + str(i) + '.arff'
  createBootstrapReplicate(trainingData, replicateFilename)
  createDecisionTree(replicateFilename)

decisionTrees_output = [""]  * numSamplings
decisionTrees_status = [True] * numSamplings
decisionTrees_running = True

while decisionTrees_running:
  for i in range(0, len(decisionTreeProcesses)):
    proc = decisionTreeProcesses[i]
    poll = proc.wait(os.WNOHANG)
    
    # if poll is not None, process is finished
    if poll is not None:
       decisionTrees_status[i] = False

    # get output of process if any
    output = proc.read()
    if output != "":
      decisionTrees_output[i] += output

  decisionTrees_running = areTreesRunning(decisionTrees_status)

decisionTrees = []
for treeString in decisionTrees_output:
  decisionTrees.append(json.loads(treeString))


# We create 20 decision trees once, and then run tests against subsets of the 20 for the smaller ensembles
runTestWithDecisionTrees(decisionTrees[0:1])
runTestWithDecisionTrees(decisionTrees[0:3])
runTestWithDecisionTrees(decisionTrees[0:5])
runTestWithDecisionTrees(decisionTrees[0:10])
runTestWithDecisionTrees(decisionTrees)
print 'Complete'