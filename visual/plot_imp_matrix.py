#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import seaborn as sns
import os
import sparsebinarydistance.distance as distance
import traceback

matplotlib.rcParams['figure.dpi']= 300


def select_subset(snvmatrix, n=0):
    snvmatrix = snvmatrix.loc[(snvmatrix >= 0).sum(axis=1) >= n]
    return snvmatrix

def gen_vis_matrix(raw_matrix, imputed_matrix, transparency = 0.45):
    '''
    raw values -> 0,1
    imputed values -> 0.45, 0.55
    return post_impute_matrix
    return a matrix with graded values representing the evidence of values (can be even more! different round of imputation!)
    '''
    # check if values only contains 0 and 1.

    print("raw == 1: " ,(raw_matrix == 1).sum().sum())
    print("raw == 0: " ,(raw_matrix == 0).sum().sum())


    # 2. enhance the transparency of the imputed values
    transparent_matrix = imputed_matrix.copy()
    transparent_matrix[(transparent_matrix == 0)] = transparency
    transparent_matrix[(transparent_matrix == 1)] = ( 1 -transparency)

    print("imputed == 0.55: ", (transparent_matrix == 1- transparency).sum().sum())
    print("imputed == 0.45: ", (transparent_matrix == transparency).sum().sum())

    # 3. combine the measured values and imputed values in a graded matrix
    # Select index from the imputed cells from the raw matrix
    raw_matrix = raw_matrix.loc[list(imputed_matrix.index)]
    print('shape of original and imputed matrix are the same:', raw_matrix.shape == transparent_matrix.shape)
    ### if the shape is not the same, doesn't update.

    combined_matrix = transparent_matrix.copy()
    combined_matrix.update(raw_matrix)

    # print("combined total 1 equals to original raw:",
    #       (raw_matrix == 1).sum().sum() == (combined_matrix == 1).sum().sum())
    # since siclonefit flip erronous observation, the following test is "False"?
    # print("combined 0.55 and 1 equals to imputed 1:", (transparent_matrix == 1 - transparency).sum().sum() == (
    #             (combined_matrix == 1).sum().sum() + (combined_matrix == 1 - transparency).sum().sum()))

    return combined_matrix



class CellSnvPlotMatrix():
    def __init__(self, raw_matrix, imputed_matrix, minPresence, minMeasurementsPerCell,
                 outpath, filename_prefix, copyNumberClone = None, replicate = None, transparency = 0.3, svg=False):
        self.raw_matrix = raw_matrix
        self.imputed_matrix = imputed_matrix
        self.sanity_input_matrix()
        self.transparency = transparency
        self.svg = svg
        self.minPresence = minPresence
        self.minMeasurementsPerCell = minMeasurementsPerCell
        self.replicate = replicate
        self.copyNumberClone = copyNumberClone

        self.clusterNumberGroup_mapping = None
        self.clusterNumber = None
        self.cnv_group_name = None
        self.gen_CNVGroup() # create copyNumberClone, clusterNumberGroup_mapping

        self.outpath = outpath
        self.filename_prefix = filename_prefix
        # form by createClusterGroupsColor() <- called in gen_mapping_color()
        self.clusterGroupColors = None
        # form by cal_sparse_distance()
        self.keptsSNVs = None
        self.keptsSNVs_for_plotting = None
        self.jointMatrixCluster = None
        self.snv_distance_matrix = None
        # form by plot_snv_snv() <- after cal_sparse_distance()
        self.snv_order = None
        # form by gen_mapping_color() <- after plot_snv_snv()
        self.jm = None
        self.colors = None
        self.mapping = None
        self.cell_order = None
        self.cal_sparse_distance()
        self.plot_snv_snv()
        self.gen_mapping_color()
        self.plot_cell_cell()


    def gen_CNVGroup(self):
        # 0 :Chr18- (Complete loss, no Chr4 loss)
        # 1: Chr18- + Chr4A -,
        # 2: Chr18-, Chr 4B-
        # 3: Chr8- (No Chr18- no chr4-)
        # 4: Chr9+ (No Chr18- no chr4- no chr8-)
        # 4: Others
        # 2020/04/13: new group total 52 states
        clusterNumberGroup = {
            0: [2, 9, 15, 17, 28, 29, 42, 47, 52],  # chr18- (Whole chr18 loss)
            1: [3, 13, 27, 34, 40, 41, 43],  # chr18- chr4A
            2: [4, 22, 36, 46],  # chr18- chr4B
            3: [1, 8, 16, 18, 19, 23, 30, 32, 35, 44, 49, 50, 51],  # chr8- (or chr8+)
            4: [5, 6, 7, 10, 12, 20, 25, 26, 31, 38, 39, 45],  # chr9+ (also don't have chr8+ or chr18- nor chr4-)
            5: [11, 14, 21, 24, 33, 37, 48]  # others
            }
        # create mapping
        clusterNumberGroup_mapping = {}
        for k, v in clusterNumberGroup.items():
            for i in v:
                clusterNumberGroup_mapping[i] = k
        self.copyNumberClone['CNV Group'] = self.copyNumberClone['state'].map(clusterNumberGroup_mapping)
        # get only two columns I need and rewrite!
        self.copyNumberClone = self.copyNumberClone[['CNV Group','state']]
        self.clusterNumberGroup_mapping = clusterNumberGroup_mapping
        self.cnv_group_name = ["chr18(-)", "chr4A(-)", "chr4B(-)", "chr8(-)", "chr9(+)", "others"]


    def cal_sparse_distance(self, calculate_imputed_distance = False):
        if calculate_imputed_distance is False:
            cal_matrix = self.raw_matrix
        else:
            cal_matrix = self.imputed_matrix
        try:
            self.keptsSNVs, self.jointMatrixCluster, cc, dd, ee = distance.sparseDistance(
                cal_matrix,
                minPresence=self.minPresence,
                minMeasurementsPerCell=self.minMeasurementsPerCell)
            self.keptsSNVs.to_pickle(f"{self.outpath}/{self.replicate}_keptsSNVs.pickle")
            self.keptsSNVs_for_plotting = gen_vis_matrix(self.keptsSNVs, self.imputed_matrix, transparency= self.transparency)
            self.keptsSNVs_for_plotting = self.keptsSNVs_for_plotting.loc[self.keptsSNVs.index, self.keptsSNVs.columns ]
            a, self.snv_distance_matrix, c, d, e = distance.sparseDistance(
                self.keptsSNVs.T,
                minPresence=1,
                minMeasurementsPerCell=1)
        except ValueError:
            traceback.print_exc()
            print(f"{self.keptsSNVs.shape}. Module sparseDistance broke.")


    def gen_mapping_color(self):
        '''
        TODO: select what categories to add, these labels should not be included in the sparse distance calculation. But should be include in plot table
        :return:
        '''
        # self.keptsSNVs_for_plotting.index.name = None
        # self.copyNumberClone = self.copyNumberClone.loc[list(self.keptsSNVs_for_plotting.index)].reindex(self.keptsSNVs_for_plotting.index)
        # pd.to_pickle(self.keptsSNVs_for_plotting, "/Users/Alice/Desktop/keptsSNVs_for_plotting.pickle")
        # pd.to_pickle(self.copyNumberClone, "/Users/Alice/Desktop/copyNumberClone.pickle")

        # BUG: if sample names not the same, break. Now: left : with P, right: without P.
        self.jm = self.keptsSNVs_for_plotting.join(self.copyNumberClone)
        self.jm['replicate'] = self.jm.index.get_level_values(0)
        self.createClusterGroupsColor()
        self.colors, self.mapping = self.createRowColorDataFrame(nanColor=(1, 1, 1),
                                                                 predeterminedColorMapping=self.clusterGroupColors)
        self.plot_colorbar()


    def plot_snv_snv(self):
        '''
        :param minPresence:
        :param minMeasurementsPerCell:
        :return: figure saved
        '''
        cm_rows = sns.clustermap(self.snv_distance_matrix,
                                 method='ward',
                                 cmap='PiYG_r')        #distance is defined in snv_distance_matrix

        cm_rows.ax_heatmap.set_xticks([])
        cm_rows.ax_heatmap.set_xlabel(f"{len(self.snv_distance_matrix.columns)} sSNVs")
        cm_rows.ax_heatmap.set_ylabel(f"{len(self.snv_distance_matrix.columns)} sSNVs")

        cm_rows.savefig(f"{self.outpath}/{self.replicate}_{self.filename_prefix}_mm{self.minMeasurementsPerCell}_"
                        f"mp{self.minPresence}_fig0_snv_snv_distance.png", dpi=300)

        # retrieve snv order for later use
        self.snv_order = cm_rows.dendrogram_row.linkage


    def plot_cell_cell(self, shift_label=False):
        if type(self.snv_order) == type(None):
            # plot_snv_snv has not yet been done.
            self.plot_snv_snv()
            self.gen_mapping_color()
        self._plot_cell_cell(shift_label=shift_label)


    def _plot_cell_cell(self, shift_label=False):

        cm = sns.clustermap( self.jointMatrixCluster,
                       row_colors=self.colors,
                       row_cluster=True,
                       cmap='PiYG_r',  vmin=0, vmax=1, method='ward') #
        cm.ax_heatmap.set_xticks([])
        cm.ax_heatmap.set_yticks([])
        cm.ax_heatmap.set_xlabel(f'{self.jointMatrixCluster.shape[0]} single cells')
        cm.ax_heatmap.set_ylabel(f'{self.jointMatrixCluster.shape[0]} single cells')

        if shift_label == True:
            self.addpatch(cm)

        cm.savefig(f'{self.outpath}/{self.replicate}_{self.filename_prefix}_mm{self.minMeasurementsPerCell}_'
                   f'mp{self.minPresence}fig1_cell_cell.png', dpi=300)

        self.cell_order = cm.dendrogram_row.linkage

    def plot_snv_cell(self, sorted = True):
        if type(self.snv_order) == type(None):
            self.plot_snv_snv() # get self.snv_order
            self.plot_cell_cell() # get self.cell_order
            self.gen_mapping_color() # get self.jm, self.mapping, self.colors
        self._plot_snv_cell(sorted=sorted)



    def _plot_snv_cell(self, sorted = True):
        '''
        need self.jm, self.snv_order
        :param sorted:
        :return:
        '''

        if sorted == True:
            # sort matrix by column CNV Group and cluster
            sorted_matrix_combined = self.jm.sort_values(by=['CNV Group', 'replicate', 'state'], ascending=[True, True, True])
            # select the cells with CNV measurement.
            sorted_matrix_combined = sorted_matrix_combined[sorted_matrix_combined['CNV Group'].notnull()]
            sorted_matrix = sorted_matrix_combined.drop(columns=['CNV Group', 'replicate', 'state'])
            self._plot_sorted(sorted_matrix)
        else:
            cg = sns.clustermap(self.keptsSNVs_for_plotting.fillna(0.5).T,
                                cmap='bwr',
                                xticklabels=False,
                                yticklabels=1,
                                # col_linkage=self.cell_order,
                                # row_linkage=self.snv_order,
                                col_colors=self.colors,
                                method='ward')
            cg.ax_heatmap.set_xlabel(f"{len(self.keptsSNVs.index)} single cells")
            cg.ax_heatmap.set_ylabel(f"{len(self.keptsSNVs.columns)} sSNV positions")
            cg.savefig(f'{self.outpath}/{self.replicate}_{self.filename_prefix}_mm{self.minMeasurementsPerCell}_'
                       f'mp{self.minPresence}_fig2_cell_snv.png', dpi=300)


    def _plot_sorted(self, sorted_matrix):
        # sorted_matrix = pd.read_pickle("/Users/Alice/Desktop/sorted_matrix.pickle")
        # passed:
        # self.snv_order = pd.read_pickle("/Users/Alice/Desktop/row_order.pickle")
        # self.colors = pd.read_pickle("/Users/Alice/Desktop/colors.pickle")

        cm_patch_sorted = sns.clustermap(sorted_matrix.fillna(0.5).T,
                                         cmap='bwr',
                                         yticklabels=False,
                                         xticklabels=False,
                                         col_cluster=False,
                                         # row_linkage=self.snv_order,
                                         col_colors=self.colors,
                                         method='ward')

        cm_patch_sorted.ax_heatmap.set_xlabel(f"{len(self.keptsSNVs.index)} single cells")
        cm_patch_sorted.ax_heatmap.set_ylabel(f"{len(self.keptsSNVs.columns)} sSNV positions")
        cm_patch_sorted.savefig(f"{self.outpath}/{self.replicate}_{self.filename_prefix}__"
                                f"mm{self.minMeasurementsPerCell}_mp{self.minPresence}_fig3_cell_snv_sorted.png", dpi=300)
        if self.svg == True:
            cm_patch_sorted.savefig(f"{self.outpath}/{self.replicate}_{self.filename_prefix}_"
                                    f"mm{self.minMeasurementsPerCell}_mp{self.minPresence}_fig3_cell_snv_sorted.svg")


    def plot_colorbar(self):

        count = 1
        for name in self.mapping:
            lut = []
            lutKeys = []
            for key in sorted(list([x for x in self.mapping[name].keys() if not pd.isnull(x)])):
                lut.append(self.mapping[name][key])
                lutKeys.append(key)
            sns.palplot(sns.color_palette(lut))
            locs, labels = plt.xticks()
            # if plotting CNV group state color, plotting integer is not informative. Get cnv_group_name variable
            if type(lutKeys[0]) == int:
                lutKeys = self.cnv_group_name
            plt.xticks(locs + 0.5, lutKeys, fontsize='12')
            plt.savefig(f"{self.outpath}/{self.replicate}_{self.filename_prefix}_cmap{count}_{name}.png", bbox_inches='tight')
            count += 1


    def createClusterGroupsColor(self):

        clusterOrder = [0, 1, 2, 3, 4]
        self.clusterNumber = len(clusterOrder)
        self.clusterGroupColors = dict(zip(clusterOrder, [(1.0, 0.8667, 0.0),
                                                     (1.0, 0.6471, 0.0),
                                                     (0.6157, 0.9373, 0.2235),
                                                     (0.2588, 0.8196, 0.9569),
                                                     (0.2549, 0.4549, 0.9569)]))

    def createRowColorDataFrame(self, nanColor=(0, 0, 0), predeterminedColorMapping={}):
        # Should look like:
        # discreteStatesDataFrame = pd.DataFrame( [['A','x'],['A','y']],index=['A','B'], columns=['First', 'Second'] )
        colorMatrix = []
        luts = {}
        discreteStatesDataFrame = self.jm[['CNV Group', 'replicate']]
        for column in discreteStatesDataFrame:
            states = [x for x in discreteStatesDataFrame[column].unique() if not pd.isnull(x)]
            undeterminedColorStates = [x for x in discreteStatesDataFrame[column].unique() if
                                       not pd.isnull(x) and not x in predeterminedColorMapping]

            cols = sns.color_palette('hls', len(undeterminedColorStates))
            # lut = { i:sns.color_palette('bright').jet(x) for i,x in zip(states, np.linspace(0,1,len(states)) )}
            lut = {state: cols[i] for i, state in enumerate(undeterminedColorStates)}
            lut.update({key: value for key, value in predeterminedColorMapping.items() if key in states})
            lut[np.nan] = nanColor
            colorMatrix.append([nanColor if pd.isnull(x) else lut[x] for x in discreteStatesDataFrame[column]])
            luts[column] = lut
        discreteColorMatrix = pd.DataFrame(colorMatrix, index=discreteStatesDataFrame.columns,
                                           columns=discreteStatesDataFrame.index).transpose()
        return discreteColorMatrix, luts

    def addpatch(self, cm):

        clusterRowIndex = 0
        for i, c in enumerate(np.array(cm.row_colors[0])[cm.dendrogram_row.reordered_ind]):

            cluster = self.jm['CNV Group'].iloc[cm.dendrogram_row.reordered_ind[i]]

            # add a white rectangular to mask the original row_color
            rect = patches.Rectangle(
                (i, clusterRowIndex),
                1, 4, fill=True, facecolor='w', edgecolor='blue', lw=0)
            cm.ax_row_colors.add_patch(rect)

            if np.isnan(cluster):
                continue

            xStep = (1 / (self.clusterNumber + 1))
            startX = xStep * self.clusterNumberGroup_mapping.get(cluster, 1)
            rect = patches.Rectangle(
                (i, clusterRowIndex + startX),
                1, xStep * 4, fill=True, facecolor=c, edgecolor='k', lw=0)
            cm.ax_row_colors.add_patch(rect)


    def sanity_input_matrix(self):
        # 1. check if only include 0, 1. If -1: convert to np.nan
        self.raw_matrix[(self.raw_matrix == -1)] = np.nan
        self.raw_matrix[(self.raw_matrix > 0.5)] = 1
        self.raw_matrix[(self.raw_matrix < 0.5)] = 0
        self.imputed_matrix[(self.imputed_matrix == -1)] = np.nan
        self.imputed_matrix[(self.imputed_matrix > 0.5)] = 1
        self.imputed_matrix[(self.imputed_matrix < 0.5)] = 0


# THis method doesn't work because this only alters
def check_indexname(pickle_files):
    '''

    :param pickle_files: a list of pickle files
    :return:
    '''
    pickle_files_new = []
    for pickle_file in pickle_files:
        if pickle_file.index.names is not [None, None, None, None]:
            pickle_file.index.names = [None, None, None, None]
        pickle_files_new.append(pickle_file)
    return pickle_files_new
