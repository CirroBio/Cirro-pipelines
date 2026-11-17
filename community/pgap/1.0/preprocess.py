#!/usr/bin/env python3

CONTACT_FIELDS = [
    "last_name",
    "first_name",
    "email",
    "organization",
    "department",
    "phone",
    "street",
    "city",
    "state",
    "postal_code",
    "country"
]


def write_text(uri: str, content: str):
    """Write a text file to the S3 location given by uri."""

    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.models.s3_path import S3Path

    path = S3Path(uri)
    assert path.valid, f"Not an S3 URI: {uri}"
    boto3.client("s3").put_object(
        Bucket=path.bucket,
        Key=path.key,
        Body=content.encode()
    )


def yaml_quote(value) -> str:
    """A YAML single-quoted scalar, in which a quote is escaped by doubling it."""
    return "'" + str(value).replace("'", "''") + "'"


def make_submol(params: dict) -> str:
    """The 'submol' YAML with the organism and contact annotations PGAP expects."""

    lines = [
        "organism:",
        f"    genus_species: {yaml_quote(params['genus_species'])}",
    ]
    if params.get("strain") is not None:
        lines.append(f"    strain: {yaml_quote(params['strain'])}")

    lines.append("contact_info:")
    lines.extend(f"    {kw}: {yaml_quote(params[kw])}" for kw in CONTACT_FIELDS)

    lines.extend([
        "authors:",
        "    - author:",
        f"        last_name: {yaml_quote(params['last_name'])}",
        f"        first_name: {yaml_quote(params['first_name'])}",
    ])
    if params.get("middle_initial"):
        lines.append(f"        middle_initial: {yaml_quote(params['middle_initial'])}")

    return "\n".join(lines)


def make_sample_sheet(genome: str, submol_uri: str) -> str:
    """The sample sheet the workflow reads, pairing the genome FASTA with its YAML."""
    return f"""fasta,yaml
{genome},{submol_uri}
"""


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    submol = make_submol(ds.params)
    ds.logger.info("YAML - submol")
    ds.logger.info(submol)

    # Write the YAML alongside the sample sheet, in the dataset's config/
    # folder (mapped in process-input.json)
    submol_uri = ds.params["sample_sheet"].rsplit("/", 1)[0] + "/submol.yaml"
    write_text(submol_uri, submol)

    sample_sheet = make_sample_sheet(ds.params["genome"], submol_uri)
    ds.logger.info("Sample sheet")
    ds.logger.info(sample_sheet)
    write_text(ds.params["sample_sheet"], sample_sheet)
