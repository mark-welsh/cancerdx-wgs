import sys
import vcfpy
from collections import OrderedDict


strelka2_vcf = vcfpy.Reader.from_path(sys.argv[1])
tumor_sample = sys.argv[2]
normal_sample = sys.argv[3]
output_vcf_name = sys.argv[4]

ad_format = OrderedDict({
    'ID': 'AD',
    'Number': 'R',
    'Type': 'Integer',
    'Description': 'Allelic depths for the ref and alt alleles in the order listed'
})

af_format = OrderedDict({
    'ID': 'AF',
    'Number': 'A',
    'Type': 'Float',
    'Description': 'Allele fractions of alternate alleles in the tumor (cancerdx calc)'
})

gt_format = OrderedDict({
    'ID': 'GT',
    'Number': '1',
    'Type': 'String',
    'Description': 'Genotype (cancerdx default for QCI)'
})

af_info = OrderedDict({
    'ID': 'AF',
    'Number': 'A',
    'Type': 'Float',
    'Description': 'Allele Frequency, for each ALT allele, in the same order as listed (cancerdx calc)'
})

header = strelka2_vcf.header
# NOTE: for some reason Strelka2 comes with an AF
#       tag in the INFO field even though no AF tag is present.
#       it's removed here and replaced by the cancerdx one
new_header = vcfpy.header_without_lines(header, [('INFO', 'AF')])
new_header.add_format_line(ad_format)
new_header.add_format_line(af_format)
new_header.add_format_line(gt_format)
new_header.add_info_line(af_info)

output_vcf = vcfpy.Writer.from_path(output_vcf_name, new_header)
for record in strelka2_vcf:
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

    ad = [int(tier1_ref_counts), int(tier1_alt_counts)]
    try:
        af = float(tier1_alt_counts/(tier1_alt_counts + tier1_ref_counts))
    except ZeroDivisionError:
        af = 0.0
    
    normal_ad = [int(normal_ref_counts), int(normal_alt_counts)]
    try:
        normal_af = float(normal_alt_counts/(normal_alt_counts + normal_ref_counts))
    except ZeroDivisionError:
        normal_af = 0.0
   
    # NOTE: "override" the vcfpy.Record.add_format() method, it doesn't work
    #       the way I want it to, so I do it myself
    # https://github.com/bihealth/vcfpy/blob/v0.12.0/vcfpy/record.py#L142-L154
    record.FORMAT.append('AD')
    record.call_for_sample[tumor_sample].data.setdefault('AD', ad)
    record.call_for_sample[normal_sample].data.setdefault('AD', normal_ad)
    
    record.FORMAT.append('AF')
    record.call_for_sample[tumor_sample].data.setdefault('AF', [af])
    record.call_for_sample[normal_sample].data.setdefault('AF', [normal_af])
    
    # NOTE: Qiagen Clinical Insight (QCI) platform requires VCF v4.2, which requires
    #       a GT tag if 'FORMAT' is used. Since the actual genotype is not
    #       evaluated by QCI, default genotypes are assigned here
    record.FORMAT.append('GT')
    record.call_for_sample[tumor_sample].data.setdefault('GT', '0/1')
    record.call_for_sample[normal_sample].data.setdefault('GT', '0/0')

    output_vcf.write_record(record)

