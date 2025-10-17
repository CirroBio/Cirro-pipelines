#!/usr/bin/env python3

import boto3
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


s3 = boto3.resource('s3')


def check_if_exists(s3_uri):

    uri = S3Path(s3_uri)

    assert uri.valid, f"Not an S3 URI: {s3_uri}"

    try:
        s3.Object(uri.bucket, uri.key).load()
    except:
        return False
    return True

ds = PreprocessDataset.from_running()

# If the user has selected NTC normalization
if ds.params.get('control_normalization', False):
    print("User has selected for control normalization")

    # Check to see if the NTC list is present (right next to the library.csv)
    ntc_file = ds.params.get('library').replace('/library.csv', '/controls.txt')
    print(f"The expected NTC file is {ntc_file}")

    if check_if_exists(ntc_file):
        print("The file does exist in S3")

        # If that NTC list is present, then set the flag in the workflow which indicates that NTC normalization should be used
        ds.add_param(
            "use_control_normalization",
            True
        )

        ds.add_param(
            "control_sgrna",
            ntc_file
        )

    else:
        print("The file does NOT exist in S3")
