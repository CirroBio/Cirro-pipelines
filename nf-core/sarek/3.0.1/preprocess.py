#!/usr/bin/env python3

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Make a wide manifest
manifest = ds.wide_samplesheet(
    index=["sampleIndex", "sample", "lane"],
    columns="read",
    values="file",
    column_prefix="fastq_"
)
assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# append metadata to file paths
samples = ds.samplesheet.set_index("sample")
manifest = manifest.assign(**{k: manifest["sample"].apply(v.get) for k, v in samples.items()})
ordering = ['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2']
manifest = manifest.reindex(columns=ordering)

# if 'lane' was not annotated on the input dataset files
if manifest['lane'].apply(pd.isnull).any():
    # 'lane' must be unique
    # use sample as the key
    lane = []
    unique_samples = manifest.value_counts('sample')
    for sample, n in unique_samples.items():
        for i in range(n):
            lane_tmp = sample + '_' + str(i)
            lane.append(lane_tmp)

    # update manifest
    manifest = manifest.assign(lane = lane)

# Run sanity checks before writing to manifest.csv
expected_columns = ['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2']
manifest_columns = manifest.axes[1]
non_canonical = [i for i in manifest_columns if i not in expected_columns]

# Drop additional columns cleanly without exit
if len(non_canonical) != 0:
    manifest = manifest.drop(columns=non_canonical)

# Check that status is 0/1
if 'status' in manifest_columns:
    status = manifest.status
    for i, row in status.items():
        assert row in [0, 1], print(
            f"ERROR: The column status must contain 0 (normal) and/or 1 (tumor).\nOffending entry:\n {manifest.iloc[i].to_frame().T}")

# Check that sex is XX/XY
if 'sex' in manifest_columns:
    sex = manifest.sex
    for i, row in sex.items():
        assert row in ['XX', 'XY', 'NA'], print(
            f"ERROR: The column sex must consist of XX, XY or NA.\nOffending entry:\n {manifest.iloc[i].to_frame().T}")

# Check that 'patient' 'sample' and 'lane' are unique
metadata = manifest[["patient", "sample", "lane"]]
for i, row in metadata.iterrows():
    patient = str(row['patient'])
    sample = str(row['sample'])
    lane = str(row['lane'])
    assert patient != sample and sample != lane, print(
        f"ERROR: patient, sample and lane must be unique.\n Offending entry:\n {metadata.iloc[i].to_frame().T}")

# Write manifest
manifest.to_csv("manifest.csv", index=None)

# JSON parameters !! revert after fix
params = ds.params

tools = params.get('tools')
assert tools is not None, "ERROR: You must select a variant calling tool."

# Annotation tool is allowed to be empty, init empty list if it is
# not 100% sure of annotation tool behavior if neither are selected i.e empty
annotation_tool = params.get('annotation_tool') or []

# Combine the two
tools = ','.join(map(str, tools + annotation_tool))

ds.add_param('tools', tools, overwrite=True)

genome = params.get('genome')

# if user does not select VEP/snpEff then annotation tool param does not exist.
# script sets it as empty list, use this to toggle deleting the param to avoid error.
if len(annotation_tool) != 0:
    ds.remove_param('annotation_tool')

# construct logic for dbNSFP & SpliceAI
# note reference genome selected
# if true, construct the appropriate parameters as file paths. 
dbnsfp_param = params.get('vep_dbnsfp') # true or null 
spliceai = params.get('vep_spliceai') # true or null 
database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
            'GATK.GRCh38': ['GRCh38', 'hg38']}

# dbNSFP
if dbnsfp_param:
    dbnsfp = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
    dbnsfp_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
    ds.add_param('dbnsfp', dbnsfp, overwrite=True)
    ds.add_param('dbnsfp_tbi', dbnsfp_tbi, overwrite=True)
    ds.add_param('dbnsfp_consequence', 'ALL', overwrite=True)

# Splice AI
if spliceai:
    spliceai_snv = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
    spliceai_snv_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
    spliceai_indel = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
    spliceai_indel_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
    ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
    ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
    ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
    ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)

# PON handling

if params.get('analysis_type') == 'Somatic Variant Calling':
    if genome == 'GATK.GRCh37':
        pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
        pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
        ds.add_param('pon', pon, overwrite=True)
        ds.add_param('pon_tbi', pon_tbi, overwrite=True)
    if genome == "GATK.GRCh38":
        pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
        pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
        ds.add_param('pon', pon, overwrite=True)
        ds.add_param('pon_tbi', pon_tbi, overwrite=True)    

ds.remove_param('analysis_type')


# concept for dynamic resource usage
# tuple makes use of 'multipliers' for [WGS,WES] respectively. 
# the three letter codes are shorthand for processes in 'process-compute.config'
# i.e CPU_LOW, MEM_LOW, HRS_LOW ETC.
wes = params.get('wes')
assay = 1 if wes is not None else 0

compute = {
    "LOW": [1, 1],
    "MED": [2, 1],
    "HGH": [4, 2],
    "ALN": [8, 4],
    "REC": [4, 2],
    "MKD": [16, 8],
    "FRB": [2, 1]
}

resources = ['CPU', 'MEM', 'HRS']

final = {}

for name, vals in compute.items():
    for it in resources:
        outname = str(it + '_' + name)
        if 'CPU' in outname:
            outval = vals[assay]*4
        elif 'MEM' in outname:
            outval = vals[assay]*30
        elif 'HRS' in outname:
            outval = vals[assay]*24/2 if 'MKD' not in outname else vals[0]*24/4

        final[outname] = int(outval)

ds.remove_param('wes')

# Produces a dict with a key-value for every PLACEHOLDER variable
for from_str, to_str in final.items():
    ds.update_compute(from_str, str(to_str), 'nextflow-override.config')

# If an intervals file was not selected, use --no_intervals
if not ds.params.get("intervals"):
    ds.add_param("no_intervals", True)

# log
ds.logger.info(ds.params)
