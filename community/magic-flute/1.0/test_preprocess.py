"""Offline tests for the MAGeCK Flute preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, so this is a plain
`python -m unittest` that needs neither. The S3 client is a stand-in that answers
head_object the way boto3 does: nothing on success, a ClientError otherwise.
"""
import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("magic_flute_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class ClientError(Exception):
    def __init__(self, code):
        super().__init__(f"An error occurred ({code}) when calling the HeadObject operation")
        self.response = {"Error": {"Code": code}}


class S3Client:
    exceptions = SimpleNamespace(ClientError=ClientError)

    def __init__(self, error_code=None):
        self.error_code = error_code
        self.calls = []

    def head_object(self, Bucket, Key):
        self.calls.append((Bucket, Key))
        if self.error_code is not None:
            raise ClientError(self.error_code)
        return {}


class NtcFileFor(unittest.TestCase):

    def test_controls_sit_beside_the_library(self):
        self.assertEqual(
            preprocess.ntc_file_for("s3://refs/crispr_libraries/brunello/library.csv"),
            "s3://refs/crispr_libraries/brunello/controls.txt",
        )


class S3ObjectExists(unittest.TestCase):

    def test_present_object(self):
        client = S3Client()
        self.assertTrue(preprocess.s3_object_exists(client, "refs", "lib/controls.txt"))
        self.assertEqual(client.calls, [("refs", "lib/controls.txt")])

    def test_missing_object(self):
        self.assertFalse(preprocess.s3_object_exists(S3Client("404"), "refs", "lib/controls.txt"))

    def test_other_errors_are_raised_rather_than_read_as_absence(self):
        for code in ("403", "500"):
            with self.subTest(code=code), self.assertRaises(ClientError):
                preprocess.s3_object_exists(S3Client(code), "refs", "lib/controls.txt")


if __name__ == "__main__":
    unittest.main()
