"""Offline tests for the Build Genome Browser (JBrowse) 1.0 preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, and the functions
tested here take an `exists` callable in place of S3, so this is a plain
`python -m unittest` that needs nothing beyond the standard library.

The file params are the values of the form's `multiple` pickers: JSON arrays of
absolute S3 URIs.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("jbrowse_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

DATA = "s3://project-bucket/datasets/ds-1/data"
VIEWER = "https://viewers.cirro.bio/render-service-worker/s3/project-bucket/datasets/ds-1/data"


def ignore(_message):
    pass


class FindFastaInLog(unittest.TestCase):

    def test_fasta_uris_are_found_and_indexes_skipped(self):
        log = (
            f"fasta : {DATA}/genome.fa\n"
            f"fai   : {DATA}/genome.fa.fai\n"
            f"other : {DATA}/genome.fasta.gz.fai\n"
            f"gz    : {DATA}/hg38.fna.gz\n"
        )
        self.assertEqual(
            preprocess.find_fasta_in_log(log),
            {f"{DATA}/genome.fa", f"{DATA}/hg38.fna.gz"},
        )

    def test_colour_codes_are_stripped(self):
        log = f"\x1b[0;32mfasta: {DATA}/genome.fasta\x1b[0m\n"
        self.assertEqual(preprocess.find_fasta_in_log(log), {f"{DATA}/genome.fasta"})


class SingleFasta(unittest.TestCase):

    def test_same_fasta_in_several_datasets(self):
        found = {"a": {f"{DATA}/genome.fa"}, "b": {f"{DATA}/genome.fa"}}
        self.assertEqual(preprocess.single_fasta(found), f"{DATA}/genome.fa")

    def test_no_fasta(self):
        self.assertIsNone(preprocess.single_fasta({}))

    def test_different_fasta_fails(self):
        found = {"a": {f"{DATA}/one.fa"}, "b": {f"{DATA}/two.fa"}}
        with self.assertRaisesRegex(ValueError, "different FASTA files"):
            preprocess.single_fasta(found)


class Names(unittest.TestCase):

    def test_assembly_name_drops_gz_and_extension(self):
        self.assertEqual(preprocess.infer_assembly_name(f"{DATA}/GRCh38.fasta.gz"), "GRCh38")
        self.assertEqual(preprocess.infer_assembly_name(f"{DATA}/mm10.fa"), "mm10")

    def test_track_name_drops_last_extension(self):
        self.assertEqual(preprocess.infer_name(f"{DATA}/sample1.vcf.gz"), "sample1.vcf")


class SelectedPaths(unittest.TestCase):

    def test_blank_entries_are_dropped(self):
        params = {"bam": [f" {DATA}/a.bam ", "", " "]}
        self.assertEqual(preprocess.selected_paths(params, "bam"), [f"{DATA}/a.bam"])

    def test_missing_or_empty_param(self):
        for params in ({}, {"bam": None}, {"bam": []}):
            self.assertEqual(preprocess.selected_paths(params, "bam"), [])


class MakeAssembly(unittest.TestCase):

    def test_fasta_with_fai(self):
        assembly = preprocess.make_assembly(
            "genome", f"{DATA}/genome.fa", {f"{DATA}/genome.fa.fai"}.__contains__, ignore
        )
        self.assertEqual(assembly, {
            "name": "genome",
            "sequence": {"url": f"{VIEWER}/genome.fa"},
            "fai": {"url": f"{VIEWER}/genome.fa.fai"},
        })

    def test_fasta_without_fai(self):
        assembly = preprocess.make_assembly("genome", f"{DATA}/genome.fa", set().__contains__, ignore)
        self.assertNotIn("fai", assembly)

    def test_no_fasta(self):
        self.assertEqual(
            preprocess.make_assembly("unknown", None, set().__contains__, ignore),
            {"name": "unknown"},
        )


class MakeTracks(unittest.TestCase):

    ASSEMBLY = {
        "name": "genome",
        "sequence": {"url": f"{VIEWER}/genome.fa"},
        "fai": {"url": f"{VIEWER}/genome.fa.fai"},
    }

    def test_every_type_with_its_index(self):
        params = {
            "gff": [f"{DATA}/genes.gff3.gz"],
            "bigwig": [f"{DATA}/cov.bigWig"],
            "vcf": [f"{DATA}/calls.vcf.gz"],
            "cram": [f"{DATA}/s1.cram"],
            "bam": [f"{DATA}/s1.bam", f"{DATA}/s2.bam"],
        }
        existing = {
            f"{DATA}/s1.bam.bai",
            f"{DATA}/s1.cram.crai",
            f"{DATA}/calls.vcf.gz.csi",
            f"{DATA}/genes.gff3.gz.tbi",
        }
        tracks = preprocess.make_tracks(params, self.ASSEMBLY, existing.__contains__, ignore)
        self.assertEqual(tracks, [
            {"type": "bam", "name": "s1", "file": {"url": f"{VIEWER}/s1.bam"},
             "index": {"url": f"{VIEWER}/s1.bam.bai"}},
            # No index beside it: the generator assumes <file>.bai
            {"type": "bam", "name": "s2", "file": {"url": f"{VIEWER}/s2.bam"}},
            {"type": "cram", "name": "s1", "file": {"url": f"{VIEWER}/s1.cram"},
             "index": {"url": f"{VIEWER}/s1.cram.crai"},
             "sequence_adapter": {"fasta": self.ASSEMBLY["sequence"], "fai": self.ASSEMBLY["fai"]}},
            # .csi is used when there is no .tbi
            {"type": "vcf", "name": "calls.vcf", "file": {"url": f"{VIEWER}/calls.vcf.gz"},
             "index": {"url": f"{VIEWER}/calls.vcf.gz.csi"}},
            {"type": "bigwig", "name": "cov", "file": {"url": f"{VIEWER}/cov.bigWig"}},
            {"type": "gff", "name": "genes.gff3", "file": {"url": f"{VIEWER}/genes.gff3.gz"},
             "index": {"url": f"{VIEWER}/genes.gff3.gz.tbi"}},
        ])

    def test_cram_without_assembly_fasta_has_no_sequence_adapter(self):
        tracks = preprocess.make_tracks(
            {"cram": [f"{DATA}/s1.cram"]}, {"name": "unknown"}, set().__contains__, ignore
        )
        self.assertNotIn("sequence_adapter", tracks[0])

    def test_no_files_selected(self):
        self.assertEqual(preprocess.make_tracks({}, self.ASSEMBLY, set().__contains__, ignore), [])


if __name__ == "__main__":
    unittest.main()
