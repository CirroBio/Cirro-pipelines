"""Offline tests for the Cell Ranger ARC preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("cellranger_arc_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_type,fastq_1,fastq_2,fastq_3\r\n"
GEX = "01,gex,s3://bucket/data/01_GEX_S1_L001_R1_001.fastq.gz,s3://bucket/data/01_GEX_S1_L001_R2_001.fastq.gz,\r\n"
ATAC = (
    "01,atac,s3://bucket/data/01_ATAC_S2_L001_R1_001.fastq.gz,s3://bucket/data/01_ATAC_S2_L001_R2_001.fastq.gz,"
    "s3://bucket/data/01_ATAC_S2_L001_R3_001.fastq.gz\r\n"
)
REFERENCES = "s3://refs-bucket"


def samplesheet(*rows):
    return preprocess.read_samplesheet(io.StringIO(HEADER + "".join(rows)))


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        rows = samplesheet(GEX, ATAC)
        self.assertEqual(rows["sample"].tolist(), ["01", "01"])
        self.assertEqual(rows["fastq_3"].tolist(), ["", "s3://bucket/data/01_ATAC_S2_L001_R3_001.fastq.gz"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No samplesheet was provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class CheckLibraries(unittest.TestCase):

    def test_gex_and_atac_lanes_pass(self):
        lane_2 = [row.replace("_L001_", "_L002_") for row in (GEX, ATAC)]
        preprocess.check_libraries(samplesheet(GEX, ATAC, *lane_2))

    def test_other_library_type_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "must be gex or atac; check sample\\(s\\) 01"):
            preprocess.check_libraries(samplesheet(GEX.replace(",gex,", ",GEX,"), ATAC))

    def test_sample_without_atac_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "sample\\(s\\) 02 have no atac rows"):
            preprocess.check_libraries(samplesheet(GEX, ATAC, GEX.replace("01,", "02,", 1)))

    def test_sample_without_gex_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "sample\\(s\\) 02 have no gex rows"):
            preprocess.check_libraries(samplesheet(GEX, ATAC, ATAC.replace("01,", "02,", 1)))

    def test_atac_row_without_r3_fails_clearly(self):
        no_r3 = ATAC.replace("s3://bucket/data/01_ATAC_S2_L001_R3_001.fastq.gz", "")
        with self.assertRaisesRegex(ValueError, "ATAC rows need all three reads.*check sample\\(s\\) 01"):
            preprocess.check_libraries(samplesheet(GEX, no_r3))

    def test_gex_row_with_r3_fails_clearly(self):
        gex_r3 = GEX.replace(",\r\n", ",s3://bucket/data/01_GEX_S1_L001_R3_001.fastq.gz\r\n")
        with self.assertRaisesRegex(ValueError, "clear FASTQ \\(R3\\) on the gex rows of sample\\(s\\) 01"):
            preprocess.check_libraries(samplesheet(gex_r3, ATAC))


if __name__ == "__main__":
    unittest.main()
