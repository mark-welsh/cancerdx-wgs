import argparse, gzip, json
import pandas as pd


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def process_sv_vcf(vcf_lines_list):
    genomic_lines = []
    for line in vcf_lines_list:
        record = dict()
        csq = list()
        if not line[0].startswith('#'):
            record = {
                "chrm": line[0],
                "start": int(line[1]),
                "id": line[3],
                "ref": line[4],
                "type": line[5]
            }
            info = line[7].split(';')
            for field in info:
                if 'CSQ' in field:
                    record['CSQ'] = field.replace("CSQ=", "").split(",")
                elif '=' in field:
                    key, value = field.split('=')
                    record[key] = value
            genomic_lines.append(record)
    return genomic_lines


parser = argparse.ArgumentParser()
parser.add_argument("--vcf_file", dest="vcf_file", required=True)
parser.add_argument("--archer_list", dest="archer_list", required=True)
parser.add_argument("--sample_id", dest="sample_id", required=True)
args = parser.parse_args()


with open(args.archer_list, 'r') as archer_in:
    archer_genes = [x.strip() for x in archer_in.readlines()]

# test if VCF is gzipped then deal with it accordingly
if is_gz_file(args.vcf_file):
    with gzip.open(args.vcf_file,'r') as gzip_vcf_in:
        vcf_lines = gzip_vcf_in.readlines()
else:
    with open(args.vcf_file, 'r') as vcf_in:
        vcf_lines = vcf_in.readlines()

vcf = process_sv_vcf([x.strip().split('\t') for x in vcf_lines])
with open('{}.sv_reformated.vcf'.format(args.sample_id), 'w') as fout:
    fout.write('CHROM\tPOS\tID\tREF\tALT\tTYPE\tARCHER\tCSQ\n')
    archer_status = []
    for record in vcf:
        for csq in record["CSQ"]:
            archer_status = []
            for gene in archer_genes:
                if '{}|'.format(gene) in csq:
                    archer_status.append(gene)
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                record["chrm"], record["start"], record["id"], record["ref"], record["type"], ', '.join(archer_status), csq.replace('|', '\t')))
