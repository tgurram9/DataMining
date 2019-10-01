# +-----------------------------------------------------------+
# |                                                           |
# |    DIVISIVE ANALYSIS CLUSTERING ALGORITHM (Divisive)      |
# |                                                           |
# +-----------------------------------------------------------+


# --------------------- LIST OF MODULE IMPORTS ----------------
import copy
from Levenshtein import distance
from scipy.cluster import hierarchy
import numpy as np
import matplotlib.pyplot as plt

# The below class is a simple implementation of a queue, to be used later whil doing BFS during divisive clustering


class Queue:

    # Constructor creates a list
    def __init__(self):
        self.queue = list()

    # Adding clusterMembers to queue
    def enqueue(self, data):
        # Checking to avoid duplicate entry (not mandatory)
        if data not in self.queue:
            self.queue.insert(0, data)
            return True
        return False

    # Removing the last element from the queue
    def dequeue(self):
        if len(self.queue) > 0:
            return self.queue.pop()
        return ("Queue Empty!")

    # Getting the size of the queue
    def size(self):
        return len(self.queue)

    # printing the clusterMembers of the queue
    def printQueue(self):
        return self.queue

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

            if(tempString != ""):
                genomes.append(tempString)
                tempString = ""

            genomeTitles.append(fileContent[i])
        else:
            tempString += fileContent[i]
    genomes.append(tempString)
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


class cluster:
    # two constructors for two occasions
    def __init__(self, num):
        # clusterMembers contains the indices of the clusterMembers in the cluster
        self.clusterMembers = set()
        self.clusterMembers.add(num)
        self.divIndex = -1

    # to get length of the cluster
    def getLength(self):
        return len(self.clusterMembers)

    # function to get the centroid of the cluster
    def getCentroid(self, simMatrix):
        num = 0.0
        for i in self.clusterMembers:
            for j in self.clusterMembers:
                if(i != j):
                    num += simMatrix[i][j]

        return num / (2 * len(self.clusterMembers))

    # function to add an element by divIndex to the cluster
    def addElement(self, number):
        self.clusterMembers.add(number)

    # function to remove an element by divIndex from the cluster
    def removeElement(self, number):
        if number in self.clusterMembers:
            self.clusterMembers.remove(number)
            return number

    # function to print the clusterMembers in the cluster in the form of a list
    def printCluster(self):
        print(list(self.clusterMembers))

    # function that returns the clusterMembers in the cluster in the form of a list
    def returnclusterMembers(self):
        return list(self.clusterMembers)

# returns average value while determining the value of the cell divIndexed by cluster1 and cluster2 in the similarity Matrix


def getAvSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueclusterMembers
    # cluster2 has many clusterValueclusterMembers
    avVal = 0
    lengthTotal = len(cluster1.clusterMembers) + len(cluster2.clusterMembers)
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            avVal += similarityMatrix[i][j]
    return float(avVal / lengthTotal)

# returns Maximum value while determining the value of the cell divIndexed by cluster1 and cluster2 in the similarity Matrix


def getMaxSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueclusterMembers
    # cluster2 has many clusterValueclusterMembers
    maxVal = -1
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            if similarityMatrix[i][j] > maxVal:
                maxVal = similarityMatrix[i][j]
    return maxVal

# returns Minimum value while determining the value of the cell divIndexed by cluster1 and cluster2 in the similarity Matrix


def getMinSimValue(cluster1, cluster2, similarityMatrix):
    # cluster1 has many clusterValueclusterMembers
    # cluster2 has many clusterValueclusterMembers
    minVal = 1000000
    for i in cluster1.clusterMembers:
        for j in cluster2.clusterMembers:
            # print "%d %d" % (i, j)
            if similarityMatrix[i][j] < minVal:
                minVal = similarityMatrix[i][j]
    return minVal

# returns that element(by divIndex) in the cluster that is the 'least stable'
# or in simpler words, that which has the highest value of dissimilarity and has a high chance of making its own cluster
# note that it pops that element out from the cluster
# if its a singleton cluster, it returns -1


def getMaxUnstableElement(cluster, similarityMatrix):
    if (len(cluster.clusterMembers) == 1):
        return -1
        # return -1 if there is only one element in the cluster
    else:
        val = 0
        ans = -1
        ansdivIndex = -1
        for i in cluster.clusterMembers:
            for j in cluster.clusterMembers:
                if(i != j):
                    val += similarityMatrix[i][j]
            if(val > ans):
                ans = val
                ansdivIndex = i
        return cluster.removeElement(i)

# function to determine if a cluster is a singleton or not.
# useful in catching the border cases and the termination cases of the DIANA Algorithm.


def clusterIsSingleton(cluster):
    if(len(cluster.clusterMembers) == 1):
        return True
    else:
        return False

# (divIndex of the old cluster and the old cluster is given in the input)
# returns the average distance of an element in the cluster from every other element in the SAME cluster


def getAverageDis(cluster, divIndex, similarityMatrix):
    ans = 0.0
    if clusterIsSingleton(cluster):
        return -1
    if divIndex in cluster.clusterMembers:
        for i in cluster.clusterMembers:
            ans += similarityMatrix[divIndex][i]
        return float(ans / (len(cluster.clusterMembers) - 1))

# (divIndex of old cluster is given, and the new cluster in the input)
# returns the average distance of an element in a cluster from another cluster in


def getDiffClusterDis(cluster, divIndex, similarityMatrix):
    ans = 0.0
    if len(cluster.clusterMembers) < 1:
        return -1
    for i in cluster.clusterMembers:
        ans += similarityMatrix[divIndex][i]
    return float(ans / len(cluster.clusterMembers))

# +--------------------------------------- CODING REFERENCE NOTES -------------------------------------------+
# | getAverageDis(oldCluster,i_old, similarityMatrix) > getDiffClusterDis(newCluster,i_old,similarityMatrix) |
# | then add i from old to newCluster                                                                        |
# +----------------------------------------------------------------------------------------------------------+

# this function splits the 'oldCluster1' and returns the two new clusters that are formed, without changing
# 'oldCluster1'


def splitCluster(oldCluster1, similarityMatrix):
    # 'oldCluster' is the local copy which is modified, and returned.
    oldCluster = copy.deepcopy(oldCluster1)
    # Do not split the cluster if it is a singleton cluster
    if(clusterIsSingleton(oldCluster)):
        return
    else:
        # the most unstable element is 'popped' out from the oldCluster and is added to the 'newCluster'
        newCluster = cluster()
        newCluster.addElement(getMaxUnstableElement(oldCluster, similarityMatrix))

        while True:
            # break if the cluster is singleton (the basecase)
            if clusterIsSingleton(oldCluster):
                break
            else:
                # the 'condition variable' comes in handy when
                # we see that we dont see any clusterMembers that can be added to the newCluster
                # from the oldCluster then, condition remains false, and we break out of the loop
                condition = False
                for i in oldCluster.clusterMembers.copy():  # since the set oldCluster.clusterMembers is changing while iteration, we cant change it
                    if i not in newCluster.clusterMembers:
                        if i in oldCluster.clusterMembers:  # only iterate over those clusterMembers that are there in oldCluster but not in newCluster
                            if(clusterIsSingleton(oldCluster)):
                                condition = False
                                break
                            if (getAverageDis(oldCluster, i, similarityMatrix) > getDiffClusterDis(newCluster, i, similarityMatrix)):
                                newCluster.addElement(i)
                                oldCluster.removeElement(i)
                                condition = True
                # if atleast one thing changed,then condition is True
                # else, condition is False, and then break out
                if(condition is False):
                    break
        return oldCluster, newCluster


if __name__ == "__main__":

    # get the path and create the matrix from it.
    datasetPath = raw_input('Enter the path for your dataset (int .txt file) : ')
    dataMatrix = getInput(str(datasetPath))
    rawData = preProcess(dataMatrix)

    # rawData[0] is genomeTitles
    # rawData[1] is genomes
    # rawData[2] is length

    option = input('------------\nEnter \"1\" if you WANT height of dendrogram to be the dissimilarity\nEnter any other number if you DO NOT WANT height of dendrogram to be the dissimilarity\n------------\nNumber : ')
    # +--------------------------------------^------------------------------------------+
    # |  put option = 0 if you DO NOT WANT height of dendrogram to be the dissimilarity |
    # |  put option = 1 if you WANT height of dendrogram to be the dissimilarity        |
    # +--------------------------------------v------------------------------------------+

    # initialise similarity matrix
    similarityMatrix = createSimilarityMatrix(rawData[1])

    # c is the parent cluster that contains ALL the clusterMembers in it
    c = cluster()
    for i in range(0, rawData[2]):
        c.addElement(i)

    # +------------- CODING NOTES ------------------+
    # | Creating the List -_-                       |
    # | take the parent cluster                     |
    # | split it                                    |
    # | append the last two to the end of the list  |
    # +---------------------------------------------+

    # dendrogramMatrix contains the Zmatrix before its converted into a numpy array
    dendrogramMatrix = []
    # BFS Queue is used to iteratively traverse the children of the parent after splitting
    BFSQueue = Queue()
    # add the parent cluster to the BFS Queue
    BFSQueue.enqueue(c)
    # current divIndex contains the next divIndex that must be assiged to the cluster which is not a singleton
    # since, singleton clusters are named as the member that they contain
    currentdivIndex = ((2 * rawData[2]) - 2)
    # initialist the divIndex of the parent cluster and update the currentdivIndex
    c.divIndex = currentdivIndex
    currentdivIndex -= 1
    i = 1
    print "-----GENERATION 0 -----"
    # print the first generation of children (the parent, basically)
    c.printCluster()
    while BFSQueue.size() != 0:
        # temp is the matrix that is iteratively added to the dendrogramMatrix
        temp = [0.0, 0.0, 0.0, 0.0]
        top = BFSQueue.dequeue()
        # nextInLine[0] and nextInLine[1] contain the children of the top of the BFS queue
        nextInLine = splitCluster(top, similarityMatrix)
        print "-----GENERATION %s -----" % (str(i))
        # print the next generation of children
        nextInLine[0].printCluster()
        nextInLine[1].printCluster()
        # If the first child is not a singleton, update its divIndex and update the currdivIndex
        # add the first child into the BFS Queue
        if(len(nextInLine[0].clusterMembers) != 1):
            nextInLine[0].divIndex = currentdivIndex
            currentdivIndex -= 1
            BFSQueue.enqueue(nextInLine[0])

        else:
            # if the first child is a singleton
            nextInLine[0].divIndex = nextInLine[0].clusterMembers.copy().pop()
            # initialise the divIndex to the only value present

        # If the second child is not a singleton, update its divIndex and update the currdivIndex
        # add the secpmd child into the BFS Queue
        if(len(nextInLine[1].clusterMembers) != 1):
            nextInLine[1].divIndex = currentdivIndex
            currentdivIndex -= 1
            BFSQueue.enqueue(nextInLine[1])
        else:
            # If the second child is a singleton
            nextInLine[1].divIndex = nextInLine[1].clusterMembers.copy().pop()
            # initialise the divIndex to the only value present
        # update the values of the temp list
        temp[0] = float(nextInLine[0].divIndex)
        temp[1] = float(nextInLine[1].divIndex)
        temp[3] = float(len(top.clusterMembers))
        temp[2] = temp[3]
        # if the option is 1, then the distance in the temp list, is made to be the dissimilarity between the cluster
        # and its parent
        if option == 1:
            temp[2] = float(getMaxSimValue(nextInLine[0], top, similarityMatrix))
        # add temp to the dendrogram Matrix
        dendrogramMatrix.append(temp)
        i += 1
    # since the matrix we generate is in a top-down fashion, we need to reverse it to use it later on
    dendrogramMatrix.reverse()
    print "\n\n"

    # +------------- DENDROGRAM PLOTTING ------------+
    # convert dendrogramMatrix into a numpy array
    Z = np.array(dendrogramMatrix)
    # print Z
    plt.figure()
    dn = hierarchy.dendrogram(Z)
    plt.show()
    # plot the dendrogram
