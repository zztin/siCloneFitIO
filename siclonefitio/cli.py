import sys
import argparse
import pandas as pd
import siclonefitio.formatting
import visual.plot_imp_matrix as vp


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Covert formatting from and to siCloneFit java format, '
                                                 'and visualize imputation')
    parser.add_argument('--siclonefit_path', '-j', required=False,
                        default="../hamimzafar-siclonefit/SiCloneFiTComplete.jar", type=str)
    parser.add_argument('-s', '--snvmatrix', required=True, type=str)
    parser.add_argument('-cn', '--copyNumberClone', type=str)  # a path
    parser.add_argument('-r', '--replicate', default=None, type=str)
    parser.add_argument('-o', '--outpath', required=True, type=str)
    parser.add_argument('-n', '--filename_prefix', required=True, type=str)
    parser.add_argument('-mp', '--minPresence', required=True, type=int)
    parser.add_argument('-mm', '--minMeasurementsPerCell', required=True, type=int)
    parser.add_argument('-t', '--transparency', help="plot transparency for imputed values", type=float, default=0.3)
    parser.add_argument('-svg', help="create svg figure", action="store_true")
    #    parser.add_argument('-bc', '--barcodeClone', type=str)  # a path

    args = parser.parse_args()

    siclonefitio.formatting.sifit_formatting(args.snvmatrix,
                                             args.siclonefit_path,
                                             minPresence=args.minPresence,
                                             minMeasurementsPerCell=args.minMeasurementsPerCell,
                                             path_to_cnv=args.copyNumberClone,
                                             cnv_column_name="state")
    snvmatrix, imputed_snvmatrix = siclonefitio.formatting.convert_siclonefit_result(args.snvmatrix,
                                                                                     args.minPresence,
                                                                                     args.minMeasurementsPerCell)

    # plotting
    vp.make_dir(args.outpath)
    cnv = pd.read_pickle(args.copyNumberClone)
    print("snvmatrix shape: ", snvmatrix.shape)
    print("imputed snvmatrix shape: ", imputed_snvmatrix.shape)
    print("cnv matrix shape", cnv.shape)
    [snvmatrix, imputed_snvmatrix] = vp.check_indexname([snvmatrix, imputed_snvmatrix])

    if args.replicate is not None:
        snvmatrix = snvmatrix.loc[[f'{args.replicate}']]
        imputed_snvmatrix = imputed_snvmatrix.loc[[f'{args.replicate}']]

    csplot = vp.CellSnvPlotMatrix(snvmatrix, imputed_snvmatrix,
                                  args.minPresence,
                                  args.minMeasurementsPerCell,
                                  args.outpath,
                                  args.filename_prefix,
                                  cnv,
                                  args.replicate,
                                  args.transparency,
                                  args.svg)
    csplot.plot_snv_cell(sorted=True)
    csplot.plot_snv_cell(sorted=False)
    csplot.keptsSNVs_for_plotting.to_pickle(f"{args.outpath}/{args.filename_prefix}_keptsSNVs_for_plotting.pickle")

    return 0


if __name__ == '__main__':
    sys.exit(main())
