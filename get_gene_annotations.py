import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from kegg_gene import GeneSet, Gene



def get_args():
    parser = argparse.ArgumentParser("Running kegg-annotation\n")
    parser.add_argument("-f", "--csv_file", help="File with genes as first column/or directory of files", required=True)
    parser.add_argument("-p", '--pattern', required=False)
    parser.add_argument("-g", "--genome", help="Genome", default="eco", required=False)
    parser.add_argument("-colnum", "--colnum", required=False)
    parser.add_argument("-DE", "--de", help='Processing de file', action='store_true', required=False)
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




#___________________

# Take de file and add annotations based on one or more reference genomes
# Given de file
# ortho file
# reference genomes with abbreviations, i.e genome_abbr = {'UTI89':'eci', 'K-12':'eco', 'CFT073':'ecc', '536':'ecp', 'UMN026':'eum'}


def add_ortho_annotations(de_analysis_file, ortho_matrix_file):
    de_df = pd.read_csv(de_analysis_file, index_col=0)
    mat_df = pd.read_csv(ortho_matrix_file, index_col=0)
    return de_df.merge(mat_df, how='left', right_index=True, left_index=True).T



def get_gene_info_from_dict(outliers_dict, genome_dict):
    gene_info = {}
    for key, val in outliers_dict.items():
        gene_info[key] = ('', '', '', '')
        for ref in genome_dict.keys():
            try:
                if val[ref] is np.nan:
                    continue
            except KeyError:
                continue
            else:
                gene_id = '{}:{}'.format(genome_dict[ref], val[ref])
                gene = Gene(gene_id)

                gene.kegg_get()

                _, gene_name, description = gene.get_short_info().split('\t')

                if gene_name or description:
                    print(gene_name)
                    gene_info[key] = (ref, val[ref], gene_name, description)
                    break
    df = pd.DataFrame.from_dict(gene_info, orient='index', columns=["Genome","Gene_ID", "Name", "Description"])
    return df

def proccess_de_results(annotation_df, de_file):
    de_df = pd.read_csv(de_file, index_col=0)
    out_name = Path(de_file).stem +"_edited.csv"
    out_file = Path(de_file).parent/out_name
    fdf = de_df.merge(annotation_df, how='left', left_index=True, right_index=True)
    fdf.to_csv(out_file)


if __name__ == "__main__":

    args = get_args().parse_args()

    if args.de:
        genome_dict = {v.split(",")[0].strip(): v.split(",")[1].strip() for v in open(args.genome)}
        de_analysis_file = args.csv_file
        ortho_matrix_file = args.pattern
        df = add_ortho_annotations(de_analysis_file, ortho_matrix_file)
        adf = get_gene_info_from_dict(df, genome_dict)
        proccess_de_results(adf, de_analysis_file)

    else:
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
