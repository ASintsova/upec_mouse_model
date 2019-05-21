import argparse
import pandas as pd
from pathlib import Path
from kegg_gene import GeneSet, Gene



def get_args():
    parser = argparse.ArgumentParser("Running kegg-annotation\n")
    parser.add_argument("-f", "--csv_file", help="File with genes as first column/or directory of files", required=True)
    parser.add_argument("-p", '--pattern', required=False)
    parser.add_argument("-g", "--genome", help="Genome", default="eco", required=False)
    parser.add_argument("-colnum", "--colnum", required=False)
    parser.add_argument("-o", "--out_dir", required=True)
    return parser


def process_file(csv_file, genome, out_dir, gene_col=0):
    df = pd.read_csv(Path(csv_file), index_col=gene_col)
    gene_list = df.index
    gs = GeneSet(gene_list, genome, out_dir)
    gs.get_info_for_each_gene()
    info_df = gs.get_info_df()
    fdf = df.merge(info_df, how='left', left_index=True, right_index=True)
    new_csv = Path(out_dir)/(Path(csv_file).stem+"_edited.csv")
    fdf.to_csv(new_csv)


def get_files(pattern, fdir):
    return [f for f in Path(fdir).iterdir() if pattern in f.name]



if __name__ == "__main__":

    args = get_args().parse_args()
    if args.colnum:
        c = args.colnum
    else:
        c = 0
    if Path(args.csv_file).is_file():
        process_file(args.csv_file, args.genome, args.out_dir, c)
    elif args.pattern:
        files = get_files(args.pattern, args.csv_file)
        for file in files:
            print(file.name)
            process_file(file, args.genome, args.out_dir, c)
    else:
        print("Wrong inputs")
