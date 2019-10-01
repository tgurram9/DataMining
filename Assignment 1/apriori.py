import csv
import sys

from itertools import chain, combinations
from collections import defaultdict
from optparse import OptionParser

# Adds a few spaces for pretty printing!


def haveSpaceString(no):
    emptyString = ""
    for i in range(0, no):
        emptyString = emptyString + " "
    return emptyString

# Returns the powerset of a list, given an iterator (as a list).


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))

# Return all the possible subsets of a list, in the form of a set.


def subsets(s):
    return map(set, powerset(s))

# Returns the support value given the transaction list and the frequency of ouccurrence.


def getSupport(count, transactionList):
    sup = float(float(count) / len(transactionList))
    return sup

# Makes the Frequency table for all the occurences.


def makeFreqTable(itemSet, transactionList, freqSet, localSet, minConfidence, Support):
    for item in itemSet:
        for transaction in transactionList:
            if item.issubset(transaction):
                freqSet[item] = freqSet[item] + 1
                localSet[item] = localSet[item] + 1


# Returns that item set whose elements satisfies the condition of support being greater than the given support.
def returnFreqItems(itemSet, transactionList, Support, freqSet, minConfidence):
    _itemSet = set()
    localSet = defaultdict(int)

    makeFreqTable(itemSet, transactionList, freqSet, localSet, minConfidence, Support)

    for item, count in localSet.items():
        support = float(float(count) / len(transactionList))
        # Generates the itemSet to be added depending on whether the support is greater than the minimum support or not.
        if getSupport(count, transactionList) >= Support:
            _itemSet.add(item)

    # Generates the file for the Frequency items and updates it here.
    rule_title = "Freq_Items_sup:%f,conf:%f.txt" % (Support, minConfidence)
    rule_file = open(rule_title, "a")
    for item in _itemSet:
        my_str = "%s %d\n" % (list(item), freqSet[item])
        rule_file.write(my_str)

    rule_file.close()

    return _itemSet

# Generates the sets of length 'n'. (or length, as given in the argument)


def genNLenSet(itemSet, length):
    lenSet = set()
    for a in itemSet:
        for b in itemSet:
            if len(a.union(b)) == length:
                lenSet.add(a.union(b))
    return lenSet


def getItemSetTransList(database):
    # Generate singleElementSets
    itemSet = set()
    transactionList = list()
    for row in database:
        transactionList.append(frozenset(row))
        for item in frozenset(row):
            itemSet.add(frozenset([item]))
        transaction = frozenset(row)
    return itemSet, transactionList

# Generates the file where the confidence rules are sorted by confidence.


def getRules(items, rules, support, confidence):
    sortedByConfidence = sorted(rules, key=lambda (rule, confidence): confidence)    # for item, support in sortedBySupport:
    IdealSpaceLen = 40  # This is the limiter for pretty printing.
    assn_title = "Assn_Rules_sup:%f,conf:%f.txt" % (support, confidence)
    assn_file = open(assn_title, "a")
    i = 1
    for rule, confidence in sortedByConfidence:
        # Outputs the Rules into the file.
        my_str1 = "Rule[%d]: %s " % (i, str(rule[0]))
        lenOfSpace = IdealSpaceLen - len(my_str1)
        my_str1 = my_str1 + haveSpaceString(lenOfSpace)
        my_str1 = my_str1 + ":----->     %s , %.5f\n" % (str(rule[1]), confidence)
        assn_file.write(my_str1)
        i = i + 1
    assn_file.close()


def dataFromFile(file_name):
    # Function to read from the CSV file and yield a generator
    data = []
    file = open(file_name, 'rU')
    for row in file:
        # Remove trailing comma
        row = row.strip().rstrip(',')
        # Removes interrmittent commas
        record = row.split(',')
        data.append(record)
    return data


def aprioriAlgo(database, support, minConfidence):

    itemSet, transactionList = getItemSetTransList(database)
    freqSet = defaultdict(int)
    freqItemSetCollection = dict()
    assocRules = dict()
    # Generate Single ItemSet.
    singleItemSet = returnFreqItems(itemSet, transactionList, support, freqSet, minConfidence)

    candidateSet = singleItemSet
    itemsPerSet = 2
    # The following loop generates N itemsets. (where N > 1)
    while(len(candidateSet) != 0):
        freqItemSetCollection[itemsPerSet - 1] = candidateSet
        candidateSet = genNLenSet(candidateSet, itemsPerSet)
        freqItemSet = returnFreqItems(candidateSet, transactionList, support, freqSet, minConfidence)
        candidateSet = freqItemSet
        itemsPerSet += 1

    returnItems = []
    # In the following loop, the items are updated.
    for numItemSet, groupItems in freqItemSetCollection.items():
        for item in groupItems:
            sup = float(float(freqSet[item]) / len(transactionList))
            returnItems.extend([(tuple(item), sup)])

    # In the following loop, the Needed Rules are generated.
    generatedRules = []
    for numItemSet, groupItems in freqItemSetCollection.items()[1:]:
        # Frozensets are used to remove reccurrences.
        for item in groupItems:
            _subsets = []
            for x in subsets(item):
                _subsets.append(frozenset(x))

            for element in _subsets:
                remain = item.difference(element)

                if len(remain) > 0:
                    # The itemSupport and element Support are calculated.
                    itemSup = float(float(freqSet[item]) / len(transactionList))
                    elementSup = float(float(freqSet[element]) / len(transactionList))
                    # The confidence is then calculated.
                    itemConfidence = float(float(float(itemSup) / float(elementSup)))
                    # All the Rules whose confidence is greater than the minimum confidence are only added.
                    if itemConfidence >= minConfidence:
                        generatedRules.append(((tuple(element), tuple(remain)), itemConfidence))

    return returnItems, generatedRules


if __name__ == "__main__":
    # Read from the below file
    data = dataFromFile('groceries.csv')
    # Get the minimum required support and confidence as input from the user.
    support = float(raw_input('Enter min Support: '))
    confidence = float(raw_input('Enter min Confidence: '))
    # Runs the apriori algorithm and generates the Support file.
    items, rules = aprioriAlgo(data, support, confidence)
    # Generates the Rules file.
    getRules(items, rules, support, confidence)
