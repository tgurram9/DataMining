# +-----------------------------------------------------------+
# |                                                           |
# |                        LOF ALGORITHM                      |
# |                                                           |
# +-----------------------------------------------------------+

# --------------------- LIST OF MODULE IMPORTS ----------------
import numpy
from copy import deepcopy
numpy.set_printoptions(threshold=numpy.nan)

# A Function that returns the manhattan distance between two multi-attributed features


def manhattanDistance(a, b):
    ans = 0
    for i in range(0, len(a)):
        ans += abs(a[i] - b[i])
    return ans

# generates the distance matrix of the Input, of size [row x row]
# return the matrix as a list of lists


def generateDistanceMatrix(Input, DistanceMat, cols, rows):
    for i in range(0, rows):
        for j in range(0, i + 1):
            if (i != j):
                DistanceMat[i][j] = manhattanDistance(Input[i, :], Input[j, :])
                DistanceMat[j][i] = DistanceMat[i][j]

    return DistanceMat.tolist()

# For easy Generation of the K-Nearest Neighbours, we are pre-sorting a copy of the distance matrix row-by-row.
# We also keep a matrix for reference. Called as the ordered Matrix.
# This will have the positions of the new items in the distance matrix after their relative order has been jumbled.


def sortDistanceMatrix(DistanceMat_List):
    length = len(DistanceMat_List[0])
    temp = []
    orderedList = [range(length) for i in range(0, length)]

    # Start the sorting here
    for i in range(0, length):
        DistanceMat_ListRow, orderedListRow = (list(t) for t in zip(*sorted(zip(DistanceMat_List[i], orderedList[i]))))
        for j in range(0, length):
            DistanceMat_List[i][j] = DistanceMat_ListRow[j]
            orderedList[i][j] = orderedListRow[j]

    return orderedList

# Gets the list of K-Nearest-Neighbours of a data object. This does not include the object itself.


def getListOfKNN(k, row_number, DistanceMat_List_Sorted, orderedList):
    num = DistanceMat_List_Sorted[k]
    k_copy = k
    for i in range(k + 1, len(DistanceMat_List_Sorted[0])):
        if num == DistanceMat_List_Sorted[row_number][i]:
            k_copy += 1
        else:
            break
    # k_copy is for orderedList
    # k is for DistanceMat_List_Sorted
    myList = []
    for i in range(0, k_copy + 1):
        myList.append(orderedList[row_number][i])
    return myList[1:]

# A simple function for reachability distance
def reachabilityDistance(k, row_numberA, row_numberB, DistanceMat_List, DistanceMat_List_Sorted):
    return max(DistanceMat_List_Sorted[row_numberB][k], DistanceMat_List[row_numberA][row_numberB])

# A simple function for calculating LRD of a data object.
def LRD(row_numberA, KNNs, k, DistanceMat_List, DistanceMat_List_Sorted):
    sum = 0
    for i in range(len(KNNs[row_numberA])):
        sum += reachabilityDistance(k, row_numberA, KNNs[row_numberA][i], DistanceMat_List, DistanceMat_List_Sorted)
    return float(len(KNNs[row_numberA])) / float(sum)

# A simple function for calculating the LOF of a data object
def LOF(row_numberA, LRDs, KNNs):
    sum = 0
    for i in range(0, len(KNNs[row_numberA])):
        sum += LRDs[i]
    return float(sum) / float(len(KNNs[row_numberA]) * LRDs[row_numberA])


if __name__ == "__main__":
    k = input("Enter the Value of k : ")
    Input = numpy.loadtxt("german.data-numeric.txt")
    # The "-1" is added below because the last column in the dataset is simple the class label.
    cols = len(Input[0]) - 1
    rows = len(Input[:, 0])
    # Create an empty matrix of Zeros for Distance Matrix
    DistanceMat = numpy.zeros([rows, rows])
    DistanceMat_List_Sorted = generateDistanceMatrix(Input, DistanceMat, cols, rows)
    DistanceMat_List = deepcopy(DistanceMat_List_Sorted)
    # DistanceMat is the distance Matrix as a numpy Array
    # DistanceMat_List is the distance Matrix as a list of lists
    # DistanceMat_List_Sorted is the SORTES distance Matrix as a list of lists
    orderedList = sortDistanceMatrix(DistanceMat_List_Sorted)
    # orderedList has the ordering of rows in it
    # DistanceMat_List_Sorted is the distance Matrix sorted by row, as a list of lists
    KNNs = []
    LRDs = []
    LOFs = []
    classifications = []
    # The above four Lists stores the lists of KNNs, LRDs, LOFs and classifications respectively
    for i in range(0, rows):
        KNNs.append(getListOfKNN(k, i, DistanceMat_List_Sorted, orderedList))
    for i in range(0, rows):
        LRDs.append(LRD(i, KNNs, k, DistanceMat_List, DistanceMat_List_Sorted))
    for i in range(0, rows):
        LOFs.append(LOF(i, LRDs, KNNs))
    for i in range(0, rows):
        if LOFs[i] >= 1:
            # Then its an outlier
            classifications.append(2)
        else:
            classifications.append(1)
    # Fill up the corresponding lists for each of the data objects.
    # -----PRECISION-----
    # For calculating the precision
    Ones_Matches = 0
    Ones_and_Twos_Matches = 0
    for i in range(0, rows):
        if(Input[i][-1] == classifications[i]):
            Ones_and_Twos_Matches += 1
            if (classifications[i] == 1):
                Ones_Matches += 1
    precision = float(Ones_Matches) / float(Ones_and_Twos_Matches)

    # ---- RECALL ----
    # For calculating the recall
    Two_One_Mismatch = 0
    for i in range(0, rows):
        if(Input[i][-1] == 2 and classifications[i] == 1):
            Two_One_Mismatch += 1
    recall = float(Ones_Matches) / float(Ones_Matches + Two_One_Mismatch)

    # ---- ACCURACY ----
    # For calculating the accuracy
    accuracy = float(Ones_and_Twos_Matches) / float(rows)

    # --- Print out the Precision, Recall and Accuracy.
    print "PRECISION  IS :{0} %".format(precision * 100)
    print "RECALL  IS :{0} %".format(recall * 100)
    print "ACCURACY IS :{0} %".format(accuracy * 100)
