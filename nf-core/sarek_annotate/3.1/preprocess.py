#!/usr/bin/env python3

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
import urllib.request
import json


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
    annotation_tool = ds.params.get('annotation_tool') or []

    # Combine the two
    tools = ','.join(map(str, annotation_tool))

    ds.add_param('tools', tools, overwrite=True)
    genome = ds.params.get('genome')
    ds.logger.info(f"Reference genome: {genome!r}")

    # annotation_tool is a Cirro UI staging param — always remove it
    ds.remove_param('annotation_tool', force=True)

    # construct logic for dbNSFP & SpliceAI
    dbnsfp_param = ds.params.get('vep_dbnsfp') # true or null
    spliceai = ds.params.get('vep_spliceai') # true or null
    # Maps iGenomes genome identifiers to the (directory, ucsc) naming used in VEP reference paths
    database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
                'GATK.GRCh38': ['GRCh38', 'hg38']}

    # dbNSFP — only available for human genomes
    if dbnsfp_param:
        if genome in database:
            dbnsfp = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
            dbnsfp_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
            ds.add_param('dbnsfp', dbnsfp, overwrite=True)
            ds.add_param('dbnsfp_tbi', dbnsfp_tbi, overwrite=True)
            ds.add_param('dbnsfp_consequence', 'ALL', overwrite=True)
            ds.logger.info(f"dbNSFP: enabled for {genome}")
        else:
            ds.logger.warning(f"No dbNSFP reference genome found for {genome} -- removing dbNSFP plugin")
            ds.remove_param('vep_dbnsfp')

    # SpliceAI — only available for human genomes; GRCm38 (mouse) has no reference data
    if spliceai:
        if genome not in database:
            ds.logger.warning(f"SpliceAI: no reference data available for genome {genome!r} -- skipping SpliceAI plugin")
        else:
            spliceai_snv = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
            spliceai_snv_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
            spliceai_indel = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
            spliceai_indel_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
            ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
            ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
            ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
            ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)
            ds.logger.info(f"SpliceAI: enabled for {genome}")

    # PON handling
    if genome == 'GATK.GRCh37':
        pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
        pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
        ds.add_param('pon', pon, overwrite=True)
        ds.add_param('pon_tbi', pon_tbi, overwrite=True)
        ds.logger.info("PON: added GRCh37 panel of normals")
    if genome == "GATK.GRCh38":
        pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
        pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
        ds.add_param('pon', pon, overwrite=True)
        ds.add_param('pon_tbi', pon_tbi, overwrite=True)
        ds.logger.info("PON: added GRCh38 panel of normals")

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


# Params set by Cirro infrastructure that must not be overridden by user-supplied extra JSON.
_PROTECTED_PARAMS = frozenset({
    "input",
    "outdir",
    "igenomes_base",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
    "step",  # hardcoded to "annotate" by process-input.json
})

# sarek_annotate uses update_compute (not add_param) for resource vars, so no passthrough needed
_CIRRO_PASSTHROUGH_PARAMS = frozenset({"optical_duplicate_pixel_distance"})

_MAX_EXTRA_JSON_BYTES = 10_000


def apply_extra_json_params(ds: PreprocessDataset):
    """Parse the extra_params_json textarea and merge user-supplied params into the workflow.

    Applied before filter_params_by_schema so that unrecognized keys are automatically
    removed in the subsequent schema validation step.
    """
    extra_json_str = (ds.params.get("extra_params_json") or "").strip()
    ds.remove_param("extra_params_json", force=True)

    if not extra_json_str:
        ds.logger.info("extra_params_json: no extra parameters provided")
        return

    if len(extra_json_str) > _MAX_EXTRA_JSON_BYTES:
        ds.logger.warning(
            f"extra_params_json: payload too large ({len(extra_json_str):,} chars), skipping"
        )
        return

    ds.logger.info("extra_params_json: parsing user-supplied JSON parameters")

    try:
        extra_params = json.loads(extra_json_str)
    except json.JSONDecodeError as e:
        ds.logger.warning(f"extra_params_json: JSON parse error ({e}), skipping")
        return

    if not isinstance(extra_params, dict):
        ds.logger.warning(
            f"extra_params_json: expected a JSON object but got {type(extra_params).__name__}, skipping"
        )
        return

    ds.logger.info(f"extra_params_json: found {len(extra_params)} parameter(s)")

    applied, skipped = 0, 0
    for key, value in extra_params.items():
        if not isinstance(key, str) or not key.strip():
            ds.logger.warning(f"extra_params_json: skipping invalid key {key!r}")
            skipped += 1
            continue
        if key in _PROTECTED_PARAMS:
            ds.logger.warning(f"extra_params_json: skipping protected parameter '{key}'")
            skipped += 1
            continue
        ds.logger.info(f"extra_params_json: applying {key}={value!r}")
        ds.add_param(key, value, overwrite=True)
        applied += 1

    ds.logger.info(
        f"extra_params_json: applied {applied} parameter(s), skipped {skipped} protected/invalid"
    )


def filter_params_by_schema(ds: PreprocessDataset):
    """Remove any params not present in the nf-core/sarek nextflow_schema.json.

    Fetches the schema for the user-selected version from GitHub and removes every
    ds.params key not listed there. Cleans up Cirro-only params (workflow_version,
    extra_params_json) and drops invalid extra JSON params before they reach Nextflow.
    """
    version = ds.params.get("workflow_version", "3.8.1")

    # Always remove — this is a Cirro-only param that must not reach Nextflow
    ds.remove_param("workflow_version", force=True)

    url = f"https://raw.githubusercontent.com/nf-core/sarek/{version}/nextflow_schema.json"
    ds.logger.info(f"filter_params_by_schema: fetching schema for nf-core/sarek {version}")

    try:
        with urllib.request.urlopen(url) as response:
            schema = json.loads(response.read().decode())
    except Exception as e:
        ds.logger.warning(f"filter_params_by_schema: could not fetch schema — {e}, skipping filter")
        return

    allowed = set()
    for section in {**schema.get("$defs", {}), **schema.get("definitions", {})}.values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"filter_params_by_schema: schema defines {len(allowed):,} parameters")

    removed = [
        key for key in list(ds.params.keys()) if key not in allowed and key not in _CIRRO_PASSTHROUGH_PARAMS
    ]
    for key in removed:
        ds.remove_param(key, force=True)

    if removed:
        ds.logger.info(f"filter_params_by_schema: removed {len(removed)} unrecognized param(s): {removed}")
    ds.logger.info(f"filter_params_by_schema: {len(ds.params)} param(s) remain")


if __name__ == "__main__":

    # Get the currently running dataset
    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek_annotate preprocess — workflow_version={ds.params.get('workflow_version', '3.8.1')!r}")

    preprocess(ds)

    apply_extra_json_params(ds)
    filter_params_by_schema(ds)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
