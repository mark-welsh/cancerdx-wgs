import sys
import vcfpy
from collections import OrderedDict


strelka2_vcf = vcfpy.Reader.from_path(sys.argv[1])
tumor_sample = sys.argv[2]
normal_sample = sys.argv[3]
output_vcf_name = sys.argv[4]

af_format = OrderedDict({
    'ID': 'AF',
    'Number': 'A',
    'Type': 'Float',
    'Description': 'Allele fractions of alternate alleles in the tumor (cancerdx calc)'
})

af_info = OrderedDict({
    'ID': 'AF',
    'Number': 'A',
    'Type': 'Float',
    'Description': 'Allele Frequency, for each ALT allele, in the same order as listed (cancerdx calc)'
})

header = strelka2_vcf.header
header.add_format_line(af_format)
header.add_info_line(af_info)

output_vcf = vcfpy.Writer.from_path(output_vcf_name, header)
for record in strelka2_vcf:
    try:
        print(record)
        alt = record.ALT[0].serialize()
        ref = record.REF
        tumor = record.call_for_sample[tumor_sample].data
        normal = record.call_for_sample[normal_sample].data

        # NOTE: calculation below taken from Strelka2 docs
        # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic-variant-allele-frequencies
        if record.is_snv():
            tier1_ref_counts = float(tumor[str(ref) + 'U'][0])
            tier1_alt_counts = float(tumor[str(alt) + 'U'][0])
            normal_ref_counts = float(normal[str(ref) + 'U'][0])
            normal_alt_counts = float(normal[str(alt) + 'U'][0])
        else:
            tier1_ref_counts = float(tumor['TAR'][0])
            tier1_alt_counts = float(tumor['TIR'][0])
            normal_ref_counts = float(normal['TAR'][0])
            normal_alt_counts = float(normal['TIR'][0])
        af = float(tier1_alt_counts/(tier1_alt_counts + tier1_ref_counts))
        try:
            normal_af = float(normal_alt_counts/(normal_alt_counts + normal_ref_counts))
        except ZeroDivisionError:
            normal_af = 0.0
   
        # NOTE: "override" the vcfpy.Record.add_format() method, it doesn't work
        #       the way I want it to, so I do it myself
        # https://github.com/bihealth/vcfpy/blob/v0.12.0/vcfpy/record.py#L142-L154
        record.FORMAT.append('AF')
        record.call_for_sample[tumor_sample].data.setdefault('AF', [af])
        record.call_for_sample[normal_sample].data.setdefault('AF', [normal_af])
        output_vcf.write_record(record)
    except UnicodeDecodeError:
        print('FAILED FOR {}'.format(record))
        continue
