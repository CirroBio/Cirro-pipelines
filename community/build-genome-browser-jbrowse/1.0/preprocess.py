#!/usr/bin/env python3

import json
import re
from functools import partial
from pathlib import PurePosixPath


# Match any S3 URI whose final extension is .fa / .fasta / .fna (optionally .gz),
# not followed by another extension (e.g. excludes .fa.fai, .fasta.gz.fai).
FASTA_URI_PATTERN = re.compile(
    r"(s3://\S+\.(?:fa|fasta|fna)(?:\.gz)?)(?=\s|$)",
    re.MULTILINE,
)
_ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")

# Track type (also its form param), log label, and the index suffixes looked for
# beside each file, in order of preference
TRACK_TYPES = [
    ("bam", "BAM", [".bai"]),
    ("cram", "CRAM", [".crai"]),
    ("vcf", "VCF", [".tbi", ".csi"]),
    ("bigwig", "BigWig", []),
    ("gff", "GFF3", [".tbi"]),
]


def parse_s3_uri(s3_uri):
    """Split s3://bucket/key into (bucket, key)."""
    without_scheme = s3_uri.removeprefix("s3://")
    bucket, _, key = without_scheme.partition("/")
    return bucket, key


def make_file_ref(s3_uri):
    url = "https://viewers.cirro.bio/render-service-worker/s3/" + s3_uri.removeprefix("s3://")
    return {"url": url}


def s3_exists(s3, s3_uri):
    bucket, key = parse_s3_uri(s3_uri)
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except Exception:
        return False


def find_index(file_uri, suffixes, exists):
    """Return the first adjacent index URI that exists, or None."""
    for suffix in suffixes:
        candidate = file_uri + suffix
        if exists(candidate):
            return candidate
    return None


def selected_paths(params, key):
    """Return the list of S3 URIs selected for a file param."""
    return [p.strip() for p in params.get(key) or [] if p.strip()]


def infer_name(s3_uri):
    """Derive a display name from the filename (strip the last dot-delimited segment)."""
    return PurePosixPath(s3_uri).name.rsplit(".", 1)[0]


def read_process_log(s3, s3_base):
    """Read the execution log of the dataset at s3_base."""
    bucket, key = parse_s3_uri(s3_base.rstrip("/") + "/artifacts/process.log")
    return (
        s3.get_object(Bucket=bucket, Key=key)["Body"]
        .read()
        .decode("utf-8", errors="replace")
    )


def find_fasta_in_log(log_text):
    """Return the set of FASTA S3 URIs found in an execution log by FASTA_URI_PATTERN."""
    return {m.group(1) for m in FASTA_URI_PATTERN.finditer(_ANSI_ESCAPE.sub("", log_text))}


def single_fasta(fasta_by_dataset):
    """The one FASTA named across the input datasets' logs, or None if none is."""
    all_fasta = set().union(*fasta_by_dataset.values()) if fasta_by_dataset else set()
    if len(all_fasta) > 1:
        raise ValueError(f"Input datasets reference different FASTA files: {all_fasta}")
    return next(iter(all_fasta), None)


def infer_assembly_name(fasta_uri):
    """The FASTA's filename without its extension (and without .gz)."""
    stem = PurePosixPath(fasta_uri).name
    if stem.endswith(".gz"):
        stem = stem[:-3]
    return PurePosixPath(stem).stem


def make_assembly(assembly_name, fasta_uri, exists, log):
    assembly = {"name": assembly_name}
    if fasta_uri:
        assembly["sequence"] = make_file_ref(fasta_uri)
        log(f"Reference FASTA: {fasta_uri}")
        fai_uri = fasta_uri + ".fai"
        if exists(fai_uri):
            assembly["fai"] = make_file_ref(fai_uri)
            log(f"FAI index found: {fai_uri}")
        else:
            log(f"FAI index not found at {fai_uri} — omitting from assembly")
    return assembly


def make_tracks(params, assembly, exists, log):
    """One track per selected file, with the index found beside it."""
    tracks = []
    for track_type, label, index_suffixes in TRACK_TYPES:
        for path in selected_paths(params, track_type):
            log(f"{label}: {path}")
            track = {"type": track_type, "name": infer_name(path), "file": make_file_ref(path)}
            if index_suffixes:
                idx = find_index(path, index_suffixes, exists)
                if idx:
                    log(f"  index: {idx}")
                    track["index"] = make_file_ref(idx)
                else:
                    log(f"  index: not found, will be inferred as {path}{index_suffixes[0]}")
            if track_type == "cram":
                if "sequence" in assembly and "fai" in assembly:
                    track["sequence_adapter"] = {
                        "fasta": assembly["sequence"],
                        "fai": assembly["fai"],
                    }
                    log("  sequence_adapter: attached from assembly FASTA")
                else:
                    log("  sequence_adapter: not set (no assembly FASTA with FAI)")
            tracks.append(track)
    return tracks


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    s3 = boto3.client("s3")
    exists = partial(s3_exists, s3)

    ds = PreprocessDataset.from_running()

    ds.logger.info("Params: " + json.dumps(ds.params, default=str))

    # Resolve reference FASTA: use the form param if provided, otherwise scan logs
    fasta_uri = ds.params.get("fasta") or None
    if fasta_uri:
        ds.logger.info(f"Using user-supplied reference FASTA: {fasta_uri}")
    else:
        # The run's input datasets are the ones the selected files come from
        input_datasets = ds.metadata.get("inputs", [])
        ds.logger.info(f"Scanning execution logs of {len(input_datasets)} input dataset(s) for reference FASTA")
        fasta_by_dataset = {}
        for d in input_datasets:
            ds.logger.info(f"Scanning dataset {d['id']} at {d['s3']}")
            paths = find_fasta_in_log(read_process_log(s3, d["s3"]))
            if paths:
                fasta_by_dataset[d["id"]] = paths
                ds.logger.info(f"Dataset {d['id']}: found FASTA candidates {paths}")
            else:
                ds.logger.info(f"Dataset {d['id']}: no FASTA found in execution log")

        fasta_uri = single_fasta(fasta_by_dataset)
        if not fasta_uri:
            ds.logger.info("No reference FASTA found in execution logs")

    # Assembly
    assembly_name = ds.params.get("assembly_name")
    if assembly_name:
        ds.logger.info(f"Using user-supplied assembly name: {assembly_name}")
    elif fasta_uri:
        assembly_name = infer_assembly_name(fasta_uri)
        ds.logger.info(f"Inferred assembly name from FASTA filename: {assembly_name}")
    assembly_name = assembly_name or "unknown"

    assembly = make_assembly(assembly_name, fasta_uri, exists, ds.logger.info)

    tracks = make_tracks(ds.params, assembly, exists, ds.logger.info)
    ds.logger.info(f"Total tracks: {len(tracks)}")

    inputs_data = {"assembly": assembly, "tracks": tracks}
    ds.logger.info("inputs.json:\n" + json.dumps(inputs_data, indent=2))

    # Write to the dataset's config/ folder (mapped in process-input.json)
    bucket, key = parse_s3_uri(ds.params["inputs"])
    s3.put_object(Bucket=bucket, Key=key, Body=json.dumps(inputs_data, indent=2))
