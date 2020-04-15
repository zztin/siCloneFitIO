import pandas as pd
import argparse






def parse_args():
    parser = argparse.ArgumentParser(description = "Convert cell to snv csv file to pickle file "
                                                   "after checking index and column names.")
    parser.add_argument('-c', '--csv', type = str) # a path
    parser.add_argument('-n', '--cnv', type = str)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    # check MultiIndex row names
    pickle_file = pd.read_csv(args.csv, index_col=[0, 1, 2, 3], header=[0, 1])
    if pickle_file.index.names is not [None, None, None, None]:
        pickle_file.index.names == [None, None, None, None]
    else:
        pd.to_pickle(pickle_file, args.csv.split(".",0)+".pickle")
    try:
        pickle_file.join(args.cnv)
    except Exception:
        print(Exception)

