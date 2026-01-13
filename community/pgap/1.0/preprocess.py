#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset


def build_samplesheet(ds: PreprocessDataset):

    # Make a YAML file with the 'submol' annotations
    # as expected by PGAP
    submol = f"""organism:
    genus_species: '{ds.params['genus_species']}'
"""

    if ds.params.get("strain") is not None:
        submol = f"""{submol}
    strain: '{ds.params['strain']}'
"""

    # Add the contact_info information
    contact_info_fields = "\n".join([
        f"    {kw}: '{ds.params[kw]}'"
        for kw in [
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
    ])
    submol = "\n".join([submol, "contact_info:", contact_info_fields])

    # Add the author information
    submol = f"""{submol}
authors:
    - author:
        last_name: '{ds.params['last_name']}'
        first_name: '{ds.params['first_name']}'
"""
    middle_initial = ds.params.get("middle_initial")
    if middle_initial is not None and len(middle_initial) > 0:
        submol = f"""{submol}
        middle_initial: '{middle_initial}'
"""

    # Remove any blank lines
    submol = "\n".join([
        line for line in submol.split("\n")
        if len(line.strip()) > 0
    ])

    ds.logger.info("YAML - submol")
    ds.logger.info(submol)

    # Write out the YAML
    with open("submol.yaml", "w") as handle:
        handle.write(submol)

    # Make the samplesheet
    sample_sheet = f"""fasta,yaml
{ds.params['genome']},submol.yaml
"""
    ds.logger.info("Sample sheet")
    ds.logger.info(sample_sheet)
    with open("sample_sheet.csv", "w") as handle:
        handle.write(sample_sheet)

    ds.add_param("sample_sheet", "sample_sheet.csv")


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    build_samplesheet(ds)
