# +-----------------------------------------------------------+
# |                                                           |
# |             AGGLOMERATIVE CLUSTERING ALGORITHM            |
# |                                                           |
# +-----------------------------------------------------------+

# --------------------- LIST OF MODULE IMPORTS ----------------
import copy
from Levenshtein import distance
from scipy.cluster import hierarchy
import numpy as np
import matplotlib.pyplot as plt

# this function opens the file specified in the filePath and breaks it into a matrix, 'fileContent' and returns it
# every row in the matrix is a new line in the file


def getInput(filePath):
    with open(filePath, "r") as inputFile:
        fileContent = inputFile.read().splitlines()
    inputFile.close()
    return fileContent

# return the class Labels, class Items and the length as three separate lists (except for length, which is an integer)
# the input is the data matrix we generated from the 'getInput()' function


def preProcess(fileContent):
    tempString = ""
    genomes = []
    genomeTitles = []
    for i in range(0, len(fileContent)):

        if fileContent[i].startswith('>'):
            # if the string starts with '>' then it is a class label (called genomeTitles)
            if(tempString != ""):
                genomes.append(tempString)
                tempString = ""

            genomeTitles.append(fileContent[i])
        else:
            # if it doesnt start with '>' it must be a genome sequence. and the the sequence is added till we hit another '>'
            tempString += fileContent[i]
    genomes.append(tempString)
    # genomeTiles == classLabels
    # genomes == class Items
    # len(genomes) == no. of class Items
    return genomeTitles, genomes, len(genomes)


# creates the similarity matrix, given the string list of genomes in 'genomes'
def createSimilarityMatrix(genomes):
    lenGen = len(genomes)
    similarityMatrix = []
    tempList = []
    for i in range(0, lenGen):
        # start filling from the ith element. the ith element is 0.
        # the i+1th element till lenGen is getSimilarity(genomes[i],genomes[j])
        # Basically, fills the top-right part of the similarity matrix.
        for j in range(0, i + 1):
            tempList.append(0)
        for j in range(i + 1, lenGen):
            # the metric we use to mention similarity (in this case, 'dissimilarity' is called 'Levenshtein distance' and we get
            # it from the distance() function)
            tempList.append(distance(genomes[i], genomes[j]))

        similarityMatrix.append(tempList)
        tempList = []
    # fill the bottom-left triangle of the similarity matrix
    for i in range(0, lenGen):
        for j in range(0, i):
            if(0 <= i):
                similarityMatrix[i][j] = similarityMatrix[j][i]
    return similarityMatrix

# prints any 2-D list of lists
# used for debugging purposes


def printSimMatrix(simMatrix):
    print "\n---- Working Matrix ----"
    leng = len(simMatrix[0])
    for i in range(0, leng):
        print simMatrix[i]
    print "------------------------\n"

# the class definition of a cluster


class cluster:
    "Defines a cluster"

    def __init__(self, index):
        # instance variables
        self.index = index  # the index that this cluster points to in the working matrix
        self.used = False
        # this is used to define if we have processed the cluster in the current generation
        self.clusterMembers = set()
        self.clusterLength = 0
        self.ZIndex = index  # the coordinate for the Zmatrix
        self.iteration = 1  # the last coordinate for the iteration

    def addToCluster(self, num):  # adds a new member to the cluster
        self.clusterMembers.add(num)
        self.clusterLength += 1

    def printCluster(self):  # prints the members of the clusuter
        print list(self.clusterMembers)

    def returnListItems(self):  # returns the members of the cluster as a list
        return list(self.clusterMembers)

# returns average value while determining the value of the cell indexed by cluster1 and cluster2 in the similarity Matrix


def getAvSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueElements
    # cluster2 has many clusterValueElements
    avVal = 0
    lengthTotal = len(cluster1.clusterMembers) + len(cluster2.clusterMembers)
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            avVal += similarityMatrix[i][j]
    return float(avVal / lengthTotal)

# returns Maximum value while determining the value of the cell indexed by cluster1 and cluster2 in the similarity Matrix


def getMaxSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueElements
    # cluster2 has many clusterValueElements
    maxVal = -1
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            if similarityMatrix[i][j] > maxVal:
                maxVal = similarityMatrix[i][j]
    return maxVal

# returns Minimum value while determining the value of the cell indexed by cluster1 and cluster2 in the similarity Matrix


def getMinSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueElements
    # cluster2 has many clusterValueElements
    minVal = 1000000
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            if similarityMatrix[i][j] < minVal:
                minVal = similarityMatrix[i][j]
    return minVal

# In the below definition, assume activeClusters are in Order.

# num == 0 for MIN
# num == 1 for MAX
# num == any other number for AVERAGE

# adds two clusters cluster1 and cluster2 and returns their union While modifying the working matrix
# the working matrix is the similarity matrix with a few of its rows and colums fused together
# activeClusters are those clusters present in the working matrix and are arranged in the list with respect to their indices in the
# working Matrix
# currRunningIndex is the index that needs to be assigned to the new cluster that is generated


def addClusters(cluster1, cluster2, workingMatrix, similarityMatrix, activeClusters, num, currRunningIndex):

    # getting the indices of the cluster which has lesser and the greater index amongst themselves
    # the new cluster's rows will replace the one indexed by [minClusterIndex]
    # and then, maxClusterIndex will be deleted (both row and column)
    # also, min Cluster's column will be deleted
    minClusterIndex = min(cluster1.index, cluster2.index)
    maxClusterIndex = max(cluster1.index, cluster2.index)
    newCluster = cluster(minClusterIndex)

    # the below chunk of code initialises the various instance variables in the newCluster object
    cluster1.used = True
    cluster2.used = True
    newCluster.clusterMembers = cluster1.clusterMembers | cluster2.clusterMembers
    newCluster.clusterLength = len(newCluster.clusterMembers)
    newCluster.used = True
    newCluster.ZIndex = currRunningIndex
    newCluster.iteration = cluster1.iteration + cluster2.iteration  # the iteration in which the new Cluster was generated

    # +-------------------------------- CODING REFERENCE NOTES ----------------------------------+
    # | BEFORE YOU DO THE BELOW CODE, ADD THE MATRICES, AND GENERATE THE INDEX, AND THEN ADD IT  |
    # | replace the minClusterIndex row and Column FIRST. (except workingMatrix[min][max] and    |
    # | workingMatrix[max][min])                                                                 |
    # +------------------------------------------------------------------------------------------+

    # Merge the values in the newCluster row and column
    for i in range(0, len(activeClusters)):
        # activeClusters[i] is the current cluster with which things are being compared
        if(activeClusters[i].index != maxClusterIndex):
            if(num == 0):
                temp = getMinSimValue(newCluster, activeClusters[i], similarityMatrix)
            elif(num == 1):
                temp = getMaxSimValue(newCluster, activeClusters[i], similarityMatrix)
            else:
                temp = getAvSimValue(newCluster, activeClusters[i], similarityMatrix)

            workingMatrix[activeClusters[i].index][minClusterIndex] = temp
            workingMatrix[minClusterIndex][activeClusters[i].index] = temp
    # Here, the two rows in the workingMatrix have been merged together
    # now, update all the values in activeClusters, whose index is greater than the max.

    # update the indices in the activeClusters once cluster1 and cluster2 gets deleted
    for i in range(0, len(activeClusters)):
        if (activeClusters[i].index > maxClusterIndex):
            activeClusters[i].index -= 1
    workingMatrix[minClusterIndex][minClusterIndex] = 0
    # remove p4 column. that is maxClusterIndex
    for i in range(0, len(activeClusters)):
        del workingMatrix[i][maxClusterIndex]
    del workingMatrix[maxClusterIndex]
    one = 0
    # if it is a singleton Cluster, then simply initialise the index with the number in the cluster itself
    for i in range(0, len(activeClusters)):
        if activeClusters[i].index == minClusterIndex:
            one = i
            break
    del(activeClusters[one])
    for i in range(0, len(activeClusters)):
        if activeClusters[i].index == maxClusterIndex:
            one = i
            break
    del(activeClusters[one])
    # delete both the clusters
    # then, insert the newCluster in the place where activeClusters[minClusterIndex] was
    activeClusters.insert(minClusterIndex, newCluster)
    return newCluster

# returns the indices in the working Matrix whose value is the least
def getMinValWorkMat(workingMatrix):

    minVal = 100000
    i_ans = 0
    j_ans = 0

    for i in range(0, len(workingMatrix)):
        for j in range(0, i):
            if(i != j):
                if(workingMatrix[i][j] < minVal):
                    minVal = workingMatrix[i][j]
                    i_ans = i
                    j_ans = j

    return i_ans, j_ans


if __name__ == "__main__":

    # get the path and create the matrix from it.
    datasetPath = raw_input('Enter the path for your dataset (int .txt file) : ')
    dataMatrix = getInput(str(datasetPath))
    rawData = preProcess(dataMatrix)
    # rawData[0] is genomeTitles
    # rawData[1] is genomes
    # rawData[2] is length
    parameter = input('-------\nEnter \"0\" for using MIN heuristic\nEnter \"1\" for using MAX heuristic\nEnter any other number for using AVG heuristic\n-------\nNumber : ')
    # +-------------^-------------+
    # | put parameter = 0 for MIN |
    # | put parameter = 1 for MAX |
    # | put parameter = 2 for AVG |
    # +-------------v-------------+
    option = input('------------\nEnter \"1\" if you WANT height of dendrogram to be the dissimilarity\nEnter any other number if you DO NOT WANT height of dendrogram to be the dissimilarity\n------------\nNumber : ')
    # +--------------------------------------^------------------------------------------+
    # |  put option = 0 if you DO NOT WANT height of dendrogram to be the dissimilarity |
    # |  put option = 1 if you WANT height of dendrogram to be the dissimilarity        |
    # +--------------------------------------v------------------------------------------+

    # initialise similarity and working matrix
    similarityMatrix = createSimilarityMatrix(rawData[1])
    workingMatrix = copy.deepcopy(similarityMatrix)

    clusterGenerations = []
    currActiveClusters = []
    # add the individual clusters to the acyiveClusters list
    for i in range(0, rawData[2]):
        c = cluster(i)
        c.addToCluster(i)
        currActiveClusters.append(c)
    print "----"
    # print the clusters
    for i in range(0, len(currActiveClusters)):
        print(currActiveClusters[i].returnListItems())

    # Zmatrix stores the matrix that the plt library uses for plotting the dendrogram
    Zmatrix = []
    newClusterIndex = rawData[2]
    # the new cluster that is formed shall have the index of no. of item sets

    while(len(workingMatrix[0]) != 1):
        # while the no. of clusters formed is not "1", keep iterating over and over again
        i, j = getMinValWorkMat(workingMatrix)
        temp = [float(currActiveClusters[i].ZIndex), float(currActiveClusters[j].ZIndex), 0, 0]
        # temp is the list that is pushed into ZMatrix
        # since currActiveClusters[i] and [j] are modified in the loop, we simply maintain copies while getting the minSim, maxSim
        # and the AvSim values
        currActiveClusters_i = copy.deepcopy(currActiveClusters[i])
        currActiveClusters_j = copy.deepcopy(currActiveClusters[j])

        newCluster = addClusters(currActiveClusters[i], currActiveClusters[j], workingMatrix, similarityMatrix, currActiveClusters, parameter, newClusterIndex)
        temp[3] = float(newCluster.iteration)
        temp[2] = temp[3]
        # if option is 1, then make the dendrogram height to be equal to the maximum of the dissimilarity b/w the new cluster formed
        # and each of the child clusters
        if option == 1:
            if parameter == 0:
                temp[2] = max(float(getMinSimValue(currActiveClusters_i, newCluster, similarityMatrix)), float(getMinSimValue(currActiveClusters_j, newCluster, similarityMatrix)))
            elif parameter == 1:
                temp[2] = max(float(getMaxSimValue(currActiveClusters_i, newCluster, similarityMatrix)), float(getMaxSimValue(currActiveClusters_j, newCluster, similarityMatrix)))
            else:
                temp[2] = max(float(getAvSimValue(currActiveClusters_i, newCluster, similarityMatrix)), float(getAvSimValue(currActiveClusters_j, newCluster, similarityMatrix)))

        Zmatrix.append(temp)

        newClusterIndex += 1
        # the new Cluster Index that can be assigned to the next new cluster is updated
        clusterGenerations.append(currActiveClusters)

        # Below code prints the result of the each iterion during clustering on the console
        print "----"
        for i in range(0, len(currActiveClusters)):
            print(currActiveClusters[i].returnListItems())

    print "----"

    for i in range(len(Zmatrix)):
        print Zmatrix[i]

    # +------------- DENDROGRAM PLOTTING ------------+
    # Convert Zmatrix into a numpy array
    Z = np.array(Zmatrix)
    # print Z
    plt.figure()
    dn = hierarchy.dendrogram(Z)
    plt.show()
    # plot the dendrogram
