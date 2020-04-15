#!/usr/bin/env python3
import pandas as pd
import numpy as np


# retrieved from https://github.com/BuysDB/IWSS
def sparseDistance( X, minPresence=1, minMeasurementsPerCell=1, weight=True ):
    """Calculate a distance matrix based on a boolean sparse cells/observations matrix
    Parameters
    ----------
    X : pandas dataframe
        sparse boolean cell by observation matrix, contains 1s 0s and NaNs for missing data
    minPresence : int
        mimumum amount of 1s or zeros present per feature, lower sums are pruned
    minMeasurementsPerCell : int
        minimum of total observations per cell
    Returns
    -------
    keptX : pandas dataframe
        Cell by observation matrix with only the cells and features which have been kept
    jointMatrix : numpy matrix
        Cell by cell distance matrix where the Similarity and Difference have been combined
    simMatrix : numpy matrix
        Similarity matrix
    distanceMatrix: numpy matrix
        Difference matrix
    normalisationFactor : numpy vector
        Weight per feature
    """

    prev= None
    selectedCells=None
    keptX = X.copy()
    while prev is None or prev[0]!=sum(selectedCells) or prev[1]!=sum(selectedColumns):
        if selectedCells is not None:
            prev = (sum(selectedCells), sum(selectedColumns))
        selectedCells = (keptX>=0).sum(axis=1) >= minMeasurementsPerCell # Select only cells/rows with at least one measurement
        # selectedColumns = ((keptX[selectedCells]>=0).sum(axis=0)>=minPresence)
        selectedColumns = ((keptX[selectedCells]==1).sum(axis=0)>=minPresence) & ((keptX[selectedCells]==0).sum(axis=0)>0)

        #print( f'We have {sum(selectedCells)} rows and {sum(selectedColumns)} X left' )
        keptX = keptX[selectedCells].loc[:,selectedColumns]

    if keptX.shape[0]<2:
        raise ValueError('Not enough data')

    pOnes = []
    pZeros = []
    for feature in keptX.columns:
        # Weights:
        column = keptX[feature]

        pOnes.append( -np.log2( ( np.sum(column==1)/len(column) )**2 ) ) #probability of two cells both having feature
        pZeros.append( -np.log2( ( np.sum(column==0)/len(column) )**2 ) )#probability of two cells not having feature,  (and we know it)

    pOnes = np.array(pOnes)
    pZeros = np.array(pZeros)
    #print( np.sum(np.isnan(pZeros)), np.sum(np.isnan(pOnes)),  np.min(pZeros), np.min(pOnes))
    #print(keptX.shape)
    rawMatrix = keptX.values

    #return rawMatrix
    iteration = 0

    # Similarity: how much do cells look alike?
    simMatrix =  np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )
    # What is the difference between the cells?
    distanceMatrix =  np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )

    #jointMatrix = np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )

    mv = int((len(simMatrix)*(len(simMatrix)-1))/2)

    for cai in range(rawMatrix.shape[0]):
        a = rawMatrix[cai,:]

        for cbi in range(rawMatrix.shape[0]):
            b = rawMatrix[cbi,:]

            # Un normalized distance
            pairwiseUnnormalizedDistance = np.logical_and( a==1, b==0 ) * \
            (pOnes + pZeros) + \
             np.logical_and( a==0, b==1 ) * (pZeros + pOnes) # For different batches the pOnes/pZeros is batch depended

            # Normalize the distance:
            normalisationFactor = np.sum( pOnes* (a==1)) + np.sum( pZeros*(a==0)) + \
                    np.sum( pOnes* (b==1)) + np.sum( pZeros* (b==0))

            pairwiseNormalizedDistance = np.sum(pairwiseUnnormalizedDistance) / (
                    normalisationFactor )
            distanceMatrix[cai, cbi] = pairwiseNormalizedDistance

            # Similarity calculation:
            sim = np.sum( (pOnes+pOnes) * np.logical_and( a==1, b==1 )) + \
                  np.sum( (pZeros+pZeros)*np.logical_and( a==0, b==0 ))

            normalisedSim = sim/normalisationFactor
            simMatrix[cai, cbi] = normalisedSim

            #joinedDistance =  (pairwiseNormalizedDistance*0.5+ (1.0-normalisedSim )*0.5)
            #jointMatrix[cai,cbi]= joinedDistance

            if cai==cbi:
                break

    for i in range(simMatrix.shape[0]):
        for j in range(simMatrix.shape[0]):
            simMatrix[j,i]=simMatrix[i,j]
            distanceMatrix[j,i]=distanceMatrix[i,j]
            if i==j:
                break

    simMatrix = pd.DataFrame(np.clip(simMatrix,0,1))
    simMatrix.index = keptX.index
    simMatrix.columns = keptX.index

    distanceMatrix = pd.DataFrame(np.clip(distanceMatrix,0,1))
    distanceMatrix.index = keptX.index
    distanceMatrix.columns = keptX.index

    jointMatrix = pd.DataFrame( (1-simMatrix)*0.5 + distanceMatrix*0.5 )
#    jointMatrix.index = keptX.index
    #jointMatrix.columns = keptX.index


    return keptX, jointMatrix, simMatrix, distanceMatrix, normalisationFactor



def sparseDistance_no_square( X, minPresence=1, minMeasurementsPerCell=1, weight=True ):
    """Calculate a distance matrix based on a boolean sparse cells/observations matrix
    Parameters
    ----------
    X : pandas dataframe
        sparse boolean cell by observation matrix, contains 1s 0s and NaNs for missing data
    minPresence : int
        mimumum amount of 1s or zeros present per feature, lower sums are pruned
    minMeasurementsPerCell : int
        minimum of total observations per cell
    Returns
    -------
    keptX : pandas dataframe
        Cell by observation matrix with only the cells and features which have been kept
    jointMatrix : numpy matrix
        Cell by cell distance matrix where the Similarity and Difference have been combined
    simMatrix : numpy matrix
        Similarity matrix
    distanceMatrix: numpy matrix
        Difference matrix
    normalisationFactor : numpy vector
        Weight per feature
    """

    prev= None
    selectedCells=None
    keptX = X.copy()
    while prev is None or prev[0]!=sum(selectedCells) or prev[1]!=sum(selectedColumns):
        if selectedCells is not None:
            prev = (sum(selectedCells), sum(selectedColumns))
        selectedCells = (keptX>=0).sum(axis=1) >= minMeasurementsPerCell # Select only cells/rows with at least one measurement
        selectedColumns = ((keptX[selectedCells]==1).sum(axis=0)>=minPresence) & ((keptX[selectedCells]==0).sum(axis=0)>0)
        print( f'We have {sum(selectedCells)} rows and {sum(selectedColumns)} X left' )
        keptX = keptX[selectedCells].loc[:,selectedColumns]

    if keptX.shape[0]<2:
        raise ValueError('Not enough data')

    pOnes = []
    pZeros = []
    for feature in keptX.columns:
        # Weights:
        column = keptX[feature]
        if weight:
            pOnes.append( -np.log2( ( np.sum(column==1)/len(column) ) ) ) #probability of two cells both having feature
            pZeros.append( -np.log2( ( np.sum(column==0)/len(column) ) ) )#probability of two cells not having feature,  (and we know it)
        else:
            pOnes.append( -np.log2( ( 0.5 ) ) ) #probability of two cells both having feature
            pZeros.append( -np.log2( (0.5) ))

    pOnes = np.array(pOnes)
    pZeros = np.array(pZeros)
    #print( np.sum(np.isnan(pZeros)), np.sum(np.isnan(pOnes)),  np.min(pZeros), np.min(pOnes))
    #print(keptX.shape)
    rawMatrix = keptX.values

    #return rawMatrix
    iteration = 0

    # Similarity: how much do cells look alike?
    simMatrix =  np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )
    # What is the difference between the cells?
    distanceMatrix =  np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )

    #jointMatrix = np.zeros( (rawMatrix.shape[0], rawMatrix.shape[0]) )

    mv = int((len(simMatrix)*(len(simMatrix)-1))/2)

    for cai in range(rawMatrix.shape[0]):
        a = rawMatrix[cai,:]

        for cbi in range(rawMatrix.shape[0]):
            b = rawMatrix[cbi,:]

            # Un normalized distance
            pairwiseUnnormalizedDistance = np.logical_and( a==1, b==0 ) * \
            (pOnes + pZeros) + \
             np.logical_and( a==0, b==1 ) * (pZeros + pOnes) # For different batches the pOnes/pZeros is batch depended

            # Normalize the distance:
            normalisationFactor = np.sum( pOnes* (a==1)) + np.sum( pZeros*(a==0)) + \
                    np.sum( pOnes* (b==1)) + np.sum( pZeros* (b==0))

            pairwiseNormalizedDistance = np.sum(pairwiseUnnormalizedDistance) / (
                    normalisationFactor )
            distanceMatrix[cai, cbi] = pairwiseNormalizedDistance

            # Similarity calculation:
            sim = np.sum( (pOnes+pOnes) * np.logical_and( a==1, b==1 )) + \
                  np.sum( (pZeros+pZeros)*np.logical_and( a==0, b==0 ))

            normalisedSim = sim/normalisationFactor
            simMatrix[cai, cbi] = normalisedSim

            #joinedDistance =  (pairwiseNormalizedDistance*0.5+ (1.0-normalisedSim )*0.5)
            #jointMatrix[cai,cbi]= joinedDistance

            if cai==cbi:
                break

    # Copy half of the matrix to the other side
    for i in range(simMatrix.shape[0]):
        for j in range(simMatrix.shape[0]):
            simMatrix[j,i]=simMatrix[i,j]
            distanceMatrix[j,i]=distanceMatrix[i,j]
            if i==j:
                break

    simMatrix = pd.DataFrame(np.clip(simMatrix,0,1))
    simMatrix.index = keptX.index
    simMatrix.columns = keptX.index

    distanceMatrix = pd.DataFrame(np.clip(distanceMatrix,0,1))
    distanceMatrix.index = keptX.index
    distanceMatrix.columns = keptX.index

    jointMatrix = pd.DataFrame( (1-simMatrix)*0.5 + distanceMatrix*0.5 )
#    jointMatrix.index = keptX.index
    #jointMatrix.columns = keptX.index


    return keptX, jointMatrix, simMatrix, distanceMatrix, normalisationFactor