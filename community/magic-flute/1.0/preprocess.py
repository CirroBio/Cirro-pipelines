#!/usr/bin/env python3


def ntc_file_for(library: str) -> str:
    """The list of non-targeting controls, kept right next to the library.csv."""
    return library.replace('/library.csv', '/controls.txt')


def s3_object_exists(s3_client, bucket: str, key: str) -> bool:
    """Only a 404 means the object is absent. Any other error, such as access denied,
    is raised rather than read as absence, which would silently skip NTC normalization."""

    try:
        s3_client.head_object(Bucket=bucket, Key=key)
    except s3_client.exceptions.ClientError as err:
        if err.response["Error"]["Code"] == "404":
            return False
        raise
    return True


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset
    from cirro.models.s3_path import S3Path

    ds = PreprocessDataset.from_running()

    # If the user has selected NTC normalization
    if ds.params.get('control_normalization', False):
        print("User has selected for control normalization")

        # Check to see if the NTC list is present (right next to the library.csv)
        ntc_file = ntc_file_for(ds.params.get('library'))
        print(f"The expected NTC file is {ntc_file}")

        uri = S3Path(ntc_file)
        assert uri.valid, f"Not an S3 URI: {ntc_file}"

        if s3_object_exists(boto3.client('s3'), uri.bucket, uri.key):
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

    # Consumed above; not a parameter the workflow declares.
    ds.remove_param("control_normalization", force=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
