"""Offline tests for the nf-core/demultiplex 1.0 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The run folder and samplesheet params are the values of the form's pickers: absolute
S3 URIs, the folder without a trailing slash.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("demux_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

RUN_FOLDER = "s3://bucket/datasets/run-1/data/240101_A00123_0001_AHXXXXXXX"
SAMPLESHEET = "s3://bucket/datasets/sheets/data/SampleSheet.csv"
PARAMS = {"run_folder": RUN_FOLDER, "samplesheet": SAMPLESHEET}


class MakeManifest(unittest.TestCase):

    def test_one_row_per_lane(self):
        manifest = preprocess.make_manifest(
            {**PARAMS, "flowcell_id": "AHXXXXXXX", "lane": "1,2"}, "My Run"
        )
        self.assertEqual(
            manifest.to_csv(index=False),
            "id,samplesheet,lane,flowcell\n"
            f"AHXXXXXXX,{SAMPLESHEET},1,{RUN_FOLDER}\n"
            f"AHXXXXXXX,{SAMPLESHEET},2,{RUN_FOLDER}\n",
        )

    def test_no_lane_demultiplexes_every_lane(self):
        manifest = preprocess.make_manifest(PARAMS, "My Run")
        self.assertEqual(manifest["lane"].tolist(), [""])

    def test_flowcell_id_falls_back_to_the_dataset_name_without_spaces(self):
        for flowcell_id in (None, ""):
            manifest = preprocess.make_manifest({**PARAMS, "flowcell_id": flowcell_id}, "My Run 1")
            self.assertEqual(manifest["id"].tolist(), ["MyRun1"])

    def test_missing_run_folder_fails_clearly(self):
        for run_folder in (None, ""):
            with self.assertRaisesRegex(ValueError, "select the sequencing run folder"):
                preprocess.make_manifest({**PARAMS, "run_folder": run_folder}, "My Run")


if __name__ == "__main__":
    unittest.main()
