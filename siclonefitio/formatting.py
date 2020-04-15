#!/usr/bin/env python3

import pandas as pd
import argparse
import os
# from pickle to sifit txt format
# path does not contain the last  "/"



def make_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def select_subset(X, minPresence=1, minMeasurementsPerCell=1):
    """Select a subset of minimum presence of a certain SNV / minimum measurement per cell
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

    prev = None
    selectedCells = None
    keptX = X.copy()
    while prev is None or prev[0] != sum(selectedCells) or prev[1] != sum(selectedColumns):
        if selectedCells is not None:
            prev = (sum(selectedCells), sum(selectedColumns))
        selectedCells = (keptX >= 0).sum(
            axis=1) >= minMeasurementsPerCell  # Select only cells/rows with at least one measurement
        selectedColumns = ((keptX[selectedCells] == 1).sum(axis=0) >= minPresence) & (
                    (keptX[selectedCells] == 0).sum(axis=0) > 0)
        # print( f'We have {sum(selectedCells)} rows and {sum(selectedColumns)} X left' )
        keptX = keptX[selectedCells].loc[:, selectedColumns]

    if keptX.shape[0] < 2:
        raise ValueError('Not enough data')

    return keptX


def _create_sifit_table(ssnvs, path, filename):
    '''
    pandas_df converting to txt
    '''
    df = ssnvs.copy()
    df[(df == -1)] = 3
    df = df.fillna(3)
    # print('fraction of data missing: ', round((df == 3).sum().sum() / (df.shape[0] * df.shape[1]), 4))
    df = df.astype(int)
    df3 = df.T.reset_index(drop=True)
    df2 = df3.reset_index()  # create a series of number implying the sites
    # print("dataframe shape: (ssnv, cellcount)", df.shape)
    df2.to_csv( f"{path}/{filename}_siclonefit_input.txt", header=None, index=False, sep=' ')
    # index = False: don't take columns and index
    return df


def _create_txt_cell_names(df, path, filename, addTag = False, clusternumber = None):
    '''

    :param df: input pickle file with siCloneFit numbering (0,1,3)
    :param clusternumber: cnv state
    :param filename: output filename
    :return:
    '''
    cellNames = ""
    if addTag == True:
        for snv in list(df.index):
            try:
                cluster = clusternumber[snv]
            except KeyError:
                cluster = 'NaN'
            aname = tuple(map(str, snv))
            cname = '_'.join(aname)
            bname = cname + f"_c{cluster}"
            cellNames += bname + " "
    else:
        for snv in list(df.index):
            aname = tuple(map(str, snv))
            cname = '_'.join(aname)
            cellNames += cname + " "
    with open(f'{path}/{filename}_cellNames.txt', 'w') as f:
        f.write(cellNames)


def _create_txt_gene_names(df, path, filename):
    '''
    from pandas df column name to list of names in txt format.
    write file in path (see inline.)
    '''
    geneNames = ""
    namecount = 0
    for j in list(df.columns):
        aname = tuple(map(str, j))
        bname = '_'.join(aname)
        geneNames += bname + " "
        namecount += 1
    # print("Count of gene names:", namecount)
    with open(f'{path}/{filename}_geneNames.txt', 'w') as f:
        f.write(geneNames)


def _create_cmd_arg(df, path, filename, siclonefit_path, fp=0.05, fn=0.05, out_dir='./', restart="1"):
    m = df.shape[0]
    n = df.shape[1]
    missing = str(round(df.isnull().sum().sum() / (df.shape[0] * df.shape[1]), 5))
    cmd = f"java -jar {siclonefit_path} -m {m} -n {n} -fp {fp} -fn {fn} -r {restart} -df 0 " \
          f"-ipMat {path}/{filename}_siclonefit_input.txt -missing {missing} " \
          f"-cellNames {path}/{filename}_cellNames.txt " \
          f"-geneNames {path}/{filename}_geneNames.txt " \
          f"-outDir {out_dir}"

    return cmd


def sifit_formatting(path_to_ssnvs_matrix,
                     siclonefit_path,
                     out_path,
                     minPresence=1 ,
                     minMeasurementsPerCell=1,
                     path_to_cnv=None,
                     cnv_column_name = "state"):

    filename = path_to_ssnvs_matrix.rsplit("/", 1)[1]
    path = out_path
    ssnvs = pd.read_pickle(path_to_ssnvs_matrix)
    filename = str(filename.split(".")[0])+f"_mp{minPresence}_mm{minMeasurementsPerCell}"
    # select subset of ssnvmatrix with at least n measurement of 0 and 1 in the matrix
    ssnvs = select_subset(ssnvs, minPresence, minMeasurementsPerCell)
    # print("While sifit_formatting: selected ssnvs_matrix shape, ", ssnvs.shape)
    pd.to_pickle(ssnvs, f"{path}/{filename}_siclonefit_raw.pickle")
    # create txt table
    df = _create_sifit_table(ssnvs, path, filename)
    _create_txt_gene_names(df, path, filename)
    if path_to_cnv:
        copyNumberCluster = pd.read_pickle(path_to_cnv)
        clusternumber = copyNumberCluster[cnv_column_name].fillna(0).astype(int)
        _create_txt_cell_names(df, path, filename, addTag = True, clusternumber = clusternumber)
    else:
        _create_txt_cell_names(df, path, filename)

    # create command line attributes
    cmd = _create_cmd_arg(ssnvs, path, filename, siclonefit_path, fp=0.001, fn=0.0001, out_dir=path)
    # out dir has to be ./ Otherwise does not create nested

    return cmd



def missing_percentage(ssnvs_matrix):
    """
    Calculate missing percentage to detect correct folder name output from siclonefit
    :return:
    """
    missing_percent = str(ssnvs_matrix.isnull().sum().sum() /
                        (ssnvs_matrix.shape[0] * ssnvs_matrix.shape[1])).split(".")[1][0:2]
    if len(missing_percent) !=2:
        missing_percent= missing_percent+ "0"
    missing = missing_percent + "p_missing_samples"
    return missing


def convert_siclonefit_result(path_to_ssnvs_matrix, out_path, minPresence= 1, minMeasurementsPerCell = 1, path_to_txt_ssnvs_imputed=None, missing = 90):
    '''

    :param path_to_ssnvs_matrix:
    :param path_to_txt_ssnvs_imputed:
    :return:
    '''

    filename = path_to_ssnvs_matrix.rsplit("/", 1)[1]
    path = out_path
    ssnvs_matrix = pd.read_pickle(path_to_ssnvs_matrix)
    filename = str(filename.split(".")[0])+f"_mp{minPresence}_mm{minMeasurementsPerCell}"
    # select subset of ssnvmatrix with at least n measurement of 0 and 1 in the matrix, update inplace.
    ssnvs_matrix = select_subset(ssnvs_matrix, minPresence, minMeasurementsPerCell)
    # print("While back converting: selected ssnvs_matrix shape, ", ssnvs_matrix.shape)
    pd.to_pickle(ssnvs_matrix,  f"{path}/{filename}_siclonefit_raw.pickle")
    # print("ssnvs_matrix shape", ssnvs_matrix.shape)
    # round
    # missing = str(round(ssnvs_matrix.isnull().sum().sum() /
    #                     (ssnvs_matrix.shape[0] * ssnvs_matrix.shape[1]), 2)).split(".")[1]+"p_missing_samples"
    # crop
    missing = missing_percentage(ssnvs_matrix)
    print(f"[formatting.py] siCloneFit output is saved at: {path}/{missing}")
    index = ssnvs_matrix.index
    columns = ssnvs_matrix.columns

    # read in imputed txt file generated by siclonefit
    imputed_snvmatrix = pd.read_csv(f"{path}/{missing}/best/best_MAP_predicted_genotype.txt", sep=' ', names=index)
    imputed_snvmatrix = imputed_snvmatrix.T
    imputed_snvmatrix.columns = columns

    pd.to_pickle(imputed_snvmatrix, f"{path}/{filename}_siclonefit_imputed.pickle")
    print(f"[formatting.py] imputed pandas table and raw pandas table saved at \
    {path}/{filename}_siclonefit_imputed.pickle & {path}/{filename}_siclonefit_raw.pickle")
    # print(f"Cell count: {ssnvs_matrix.shape[0]}, sSNV count: {ssnvs_matrix.shape[1]}")

    return ssnvs_matrix, imputed_snvmatrix
