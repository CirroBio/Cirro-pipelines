"""Offline tests for the Spatial Data Studio preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, so this is a plain
`python -m unittest` that needs nothing else.

The picked folders arrive as one comma-separated string of S3 URIs, a whole dataset as
its data/ folder with a trailing slash.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("spatial_data_studio_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

XENIUM = {"id": "aaaaaaaa-1111", "name": "Xenium run", "dataPath": "s3://bucket/datasets/aaaaaaaa-1111/data"}
VISIUM = {"id": "bbbbbbbb-2222", "name": "Visium run", "dataPath": "s3://bucket/datasets/bbbbbbbb-2222/data/"}


class MakeInputMap(unittest.TestCase):

    def test_whole_datasets_are_named_for_the_dataset(self):
        self.assertEqual(
            preprocess.make_input_map(
                "s3://bucket/datasets/aaaaaaaa-1111/data/,s3://bucket/datasets/bbbbbbbb-2222/data/",
                [XENIUM, VISIUM],
            ),
            {
                "Xenium run": "s3://bucket/datasets/aaaaaaaa-1111/data",
                "Visium run": "s3://bucket/datasets/bbbbbbbb-2222/data",
            },
        )

    def test_subfolder_keeps_its_path_under_the_dataset_name(self):
        self.assertEqual(
            preprocess.make_input_map(
                "s3://bucket/datasets/aaaaaaaa-1111/data/slide_1/outs/",
                [XENIUM],
            ),
            {"Xenium run/slide_1/outs": "s3://bucket/datasets/aaaaaaaa-1111/data/slide_1/outs"},
        )

    def test_folder_name_sharing_the_data_prefix_is_not_matched(self):
        with self.assertRaisesRegex(ValueError, "not inside any of the run's input datasets"):
            preprocess.make_input_map("s3://bucket/datasets/aaaaaaaa-1111/data2/", [XENIUM])

    def test_repeated_dataset_names_are_suffixed_with_their_ids(self):
        twin = dict(VISIUM, name="Xenium run")
        self.assertEqual(
            preprocess.make_input_map(
                "s3://bucket/datasets/aaaaaaaa-1111/data/,s3://bucket/datasets/bbbbbbbb-2222/data/",
                [XENIUM, twin],
            ),
            {
                "Xenium run-aaaaaaaa": "s3://bucket/datasets/aaaaaaaa-1111/data",
                "Xenium run-bbbbbbbb": "s3://bucket/datasets/bbbbbbbb-2222/data",
            },
        )

    def test_folder_picked_twice_is_mapped_once(self):
        self.assertEqual(
            preprocess.make_input_map(
                "s3://bucket/datasets/aaaaaaaa-1111/data/,s3://bucket/datasets/aaaaaaaa-1111/data/",
                [XENIUM],
            ),
            {"Xenium run": "s3://bucket/datasets/aaaaaaaa-1111/data"},
        )

    def test_no_folders_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No input folders were provided"):
                preprocess.make_input_map(value, [XENIUM])

    def test_dataset_without_a_name_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "has no name"):
            preprocess.make_input_map(
                "s3://bucket/datasets/aaaaaaaa-1111/data/", [dict(XENIUM, name="")]
            )


if __name__ == "__main__":
    unittest.main()
