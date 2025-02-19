import os 
import pandas as pd
import glob

parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
input_folder = parent_dir +"/data/all_samples"
print(input_folder)

sv_files = glob.glob(os.path.join(input_folder, "*.snp.vcf"))

dfs = []

for sv_file in sv_files:
    sample_name = os.path.basename(sv_file).split('.')[0]
    
    with open(sv_file, "r") as f:
        header_lines = [line for line in f if line.startswith("#")]

    df = pd.read_csv(sv_file, comment = '#', sep = '\t', header = None)
    
    no_of_columns = 10
    if df.shape[1] >= no_of_columns:
        df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GENOTYPE-INFO"] + [F"EXTRA{i}" for i in range (df.shape[1] - no_of_columns)]
    
    df["sample"] = sample_name

    df = df.drop(columns = ["CHROM"])
    cols = df.columns.tolist()
    cols.remove("sample")
    cols.insert(1, "sample")
    df = df[cols]

    dfs.append(df)

    final_df = pd.concat(dfs, ignore_index = True)

    output_file = input_folder + "/concatenated-snp.csv"

    with open(output_file, "w") as out_file:
        final_df.to_csv(output_file, sep = "\t", index = False)

print("File concatenated")  