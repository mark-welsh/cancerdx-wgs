import argparse, json
import pandas as pd


def process_cnv_vcf(vcf_lines_list):
    genomic_lines = []
    for line in vcf_lines_list:
        record = dict()
        if not line[0].startswith('#'):
            record = {
                "chrm": line[0],
                "start": int(line[1]),
            }
            info = line[7].split(';')
            for field in info:
                if '=' in field:
                    key, value = field.split('=')
                    record[key] = value
            genomic_lines.append(record)
    return genomic_lines


parser = argparse.ArgumentParser()
parser.add_argument("--cnr_file", dest="cnr_file", required=True)
parser.add_argument("--vcf_file", dest="vcf_file", required=True)
parser.add_argument("--sample_id", dest="sample_id", required=True)
args = parser.parse_args()

cnr_df = pd.read_csv(args.cnr_file, sep='\t')

with open(args.vcf_file, 'r') as vcf_in:
    vcf_lines = vcf_in.readlines()

vcf_df = pd.DataFrame(process_cnv_vcf([x.strip().split('\t') for x in vcf_lines]))
merged = vcf_df.merge(cnr_df, how='left', left_on=["chrm", "start"], right_on=["chromosome", "start"], indicator=True)
merged.drop(columns=["CIPOS", "CIEND", "chromosome", "end", "log2", "_merge"], inplace=True)
merged.to_csv('{}.cnv_reformated.tsv'.format(args.sample_id), sep='\t', index=False)
