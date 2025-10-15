#!/usr/bin/env python3

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

def preprocess(ds):

    # Make a manifest file as expected by nf-core/sarek:
    # Columns: patient,sample,vcf

    # The logic we use for processing is:
    #   1. Analyze any file which ends with .vcf;
    #   2. The sample name will be used for 'patient';
    #   3. If the filename starts with the sample name, remove it. Use any other text before .vcf.gz as 'sample';
    #   4. If the filename does not start with the sample name, used everything before .vcf.gz as 'sample'.

    def format_manifest(r: pd.Series):

        for cname in ['sample', 'file']:
            assert cname in r.index.values, f"Column missing: {cname}"

        # Set up the values which will be returned
        dat = dict(patient=None, sample=None, vcf=r['file'])

        # Get the filename
        fn = r['file'].rsplit("/", 1)[-1]

        # If the file does not end with .vcf.gz
        if not fn.endswith('.vcf.gz'):

            # Return an empty value
            # Rows with empty values will be filtered out in the next step
            pass

        # The file does indeed end with .vcf.gz
        else:

            # Make sure that there is a valid sample value
            assert not pd.isnull(r['sample']), f"Null value for sample ({r.to_dict()})"

            # Assign 'patient' with the value from 'sample' (removing any parent folders)
            dat['patient'] = str(r['sample']).rsplit("/", 1)[-1]

            # Make a sample name which is based on the file name (without the .vcf.gz)
            sample_name = fn[:-len(".vcf.gz")]

            # If the sample name and the patient name are the same
            if dat['patient'] == sample_name:

                # Use the same sample name as patient name
                dat['sample'] = sample_name

            # Otherwise, if the sample name starts with the patient name and has additional text
            elif sample_name.startswith(dat['patient']):

                # Make a sample name which does not include the patient ID
                sample_name = sample_name[len(dat['patient']):]

                # Remove any leading '.' or '_'
                while sample_name.startswith(('.', '_')):
                    sample_name = sample_name[1:]

                # If there is no text left
                if len(sample_name) == 0:

                    # Use the patient name
                    dat['sample'] = dat['patient']

                # But if there is indeed text remaining
                else:

                    # Use that text as the sample name
                    dat['sample'] = sample_name

            # Lastly, if the file name does not start with the sample name
            else:

                # Just use the file name as the 'sample'
                dat['sample'] = sample_name

        return pd.Series(dat)

    manifest = ds.files.apply(format_manifest, axis=1)
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # Filter out any rows with null values
    filtered_manifest = manifest.dropna()
    assert filtered_manifest.shape[0] > 0, "No files found with the .vcf.gz extension"

    print(f"{filtered_manifest.shape[0]:,} / {manifest.shape[0]:,} files passed the filter for *.vcf.gz")

    # Write out the manifest
    filtered_manifest.reindex(
        columns=["patient", "sample", "vcf"]
    ).to_csv(
        "manifest.csv",
        index=None
    )

    # Annotation tool is allowed to be empty, init empty list if it is
    # not 100% sure of annotation tool behavior if neither are selected i.e empty
    annotation_tool = ds.params.get('annotation_tool') or []

    # Combine the two
    tools = ','.join(map(str, annotation_tool))

    ds.add_param('tools', tools, overwrite=True)
    genome = ds.params.get('genome')

    # if user does not select VEP/snpEff then annotation tool param does not exist.
    # script sets it as empty list, use this to toggle deleting the param to avoid error.
    if len(annotation_tool) != 0:
        ds.remove_param('annotation_tool')

    # construct logic for dbNSFP & SpliceAI
    # note reference genome selected
    # if true, construct the appropriate parameters as file paths. 
    dbnsfp_param = ds.params.get('vep_dbnsfp') # true or null 
    spliceai = ds.params.get('vep_spliceai') # true or null 
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


    # Use the same approach for resource allocation, but simply use the lower-level value set up for WES
    # concept for dynamic resource usage
    # tuple makes use of 'multipliers' for [WGS,WES] respectively. 
    # the three letter codes are shorthand for processes in 'process-compute.config'
    # i.e CPU_LOW, MEM_LOW, HRS_LOW ETC.
    assay = 1

    compute = {"LOW": [1,1],
            "MED": [2,1],
            "HGH": [4,2],
            "ALN": [8,4],
            "REC": [4,2],
            "MKD": [16,8],
            "FRB": [2,1]}

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
            
    ## Produces a dict with a key-value for every PLACEHOLDER variable in process-compute.config.
    for from_str, to_str in final.items():
        ds.update_compute(from_str, str(to_str), 'nextflow-override.config')

    # log
    ds.logger.info(ds.params)

if __name__ == "__main__":

    # Get the currently running dataset
    ds = PreprocessDataset.from_running()

    preprocess(ds)
