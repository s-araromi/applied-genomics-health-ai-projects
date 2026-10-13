#!/usr/bin/env python3
"""Retrieve the pinned genuine Phase 3 locus without downloading the full VCF.

Python's standard library handles HTTPS and the Tabix/BGZF container. Statistical
analysis, reporting and tests are implemented in R. No genotypes are generated.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import struct
import urllib.request
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = "https://1000genomes.s3.amazonaws.com/release/20130502/"
VCF = BASE + "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
PANEL = BASE + "integrated_call_samples_v3.20130502.ALL.panel"
POSITION = 136608646


def fetch(url: str, byte_range: tuple[int, int] | None = None) -> bytes:
    """Bound response size before reading; a server ignoring Range is rejected."""
    headers = {"User-Agent": "applied-genomics-health-ai-projects/Project08"}
    limit = 1_000_000
    if byte_range:
        start, end = byte_range
        headers["Range"] = f"bytes={start}-{end}"
        limit = end - start + 1
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=90) as response:
        if byte_range:
            match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",
                                 response.headers.get("Content-Range", ""))
            if response.status != 206 or not match:
                raise ValueError("Server did not honor a bounded byte-range request")
            if tuple(map(int, match.groups()[:2])) != byte_range:
                raise ValueError("Returned byte range differs from request")
        data = response.read(limit + 1)
        if len(data) > limit or (byte_range and len(data) != limit):
            raise ValueError("Unexpected source response size")
        return data


def complete_bgzf_members(data: bytes) -> bytes:
    """Decode complete gzip members; ignore only an incomplete final range block."""
    decoded = []
    while data:
        block = zlib.decompressobj(16 + zlib.MAX_WBITS)
        output = block.decompress(data)
        if not block.eof:
            break
        decoded.append(output)
        data = block.unused_data
    if not decoded:
        raise ValueError("No complete BGZF blocks returned")
    return b"".join(decoded)


def locus_range(compressed_index: bytes) -> tuple[int, int]:
    """Use the standard Tabix linear index's 16,384-base genomic windows.

    Linear offsets are lower bounds, not necessarily exact record starts. Read
    through the next distinct bound plus one maximal BGZF block, then select the
    exact variant from decoded lines. This implementation targets a simple SNP.
    """
    raw = gzip.decompress(compressed_index)
    if raw[:4] != b"TBI\x01":
        raise ValueError("Expected a Tabix TBI index")
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, name_size = struct.unpack_from("<8i", raw, 4)
    if (n_ref, fmt, col_seq, col_beg, col_end, meta, skip) != (1, 2, 1, 2, 0, 35, 0):
        raise ValueError("Unexpected chromosome-specific VCF index schema")
    if raw[36:36 + name_size] != b"2\x00":
        raise ValueError("Index is not for chromosome 2")
    offset = 36 + name_size
    n_bins = struct.unpack_from("<i", raw, offset)[0]
    offset += 4
    for _ in range(n_bins):
        _, n_chunks = struct.unpack_from("<Ii", raw, offset)
        if n_chunks < 0:
            raise ValueError("Invalid Tabix chunk count")
        offset += 8 + 16 * n_chunks
    n_intervals = struct.unpack_from("<i", raw, offset)[0]
    offset += 4
    intervals = struct.unpack_from(f"<{n_intervals}Q", raw, offset)
    window = (POSITION - 1) >> 14
    start = intervals[window] >> 16
    end_bound = next((v >> 16 for v in intervals[window + 1:] if v >> 16 > start), None)
    if not start or end_bound is None or end_bound - start > 1_000_000:
        raise ValueError("Could not establish a bounded source interval")
    return start, end_bound + 65535


def retrieve() -> dict[str, bytes]:
    headers = complete_bgzf_members(fetch(VCF, (0, 65535))).splitlines(keepends=True)
    header_lines = [line for line in headers if line.startswith(b"#")]
    if sum(line.startswith(b"#CHROM\t") for line in header_lines) != 1:
        raise ValueError("Initial range did not contain the complete VCF header")
    interval = locus_range(fetch(VCF + ".tbi"))
    decoded = complete_bgzf_members(fetch(VCF, interval))
    records = [line for line in decoded.splitlines(keepends=True)
               if line.split(b"\t")[:5] == [b"2", b"136608646", b"rs4988235", b"G", b"A"]]
    if len(records) != 1 or not records[0].endswith(b"\n"):
        raise ValueError("Expected exactly one complete matching SNP record")
    return {"data/rs4988235.phase3.vcf": b"".join(header_lines) + records[0],
            "data/phase3_samples.panel": fetch(PANEL)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="Re-fetch both pinned source files and verify before replacing")
    args = parser.parse_args()
    with (ROOT / "data/source_manifest.csv").open(newline="") as source:
        manifest = list(csv.DictReader(source))
    expected = {row["file"]: row for row in manifest}
    if not args.refresh and all((ROOT / name).is_file() for name in expected):
        content = {name: (ROOT / name).read_bytes() for name in expected}
        action = "Verified cached source"
    else:
        content = retrieve()
        action = "Retrieved and verified source"
    # Validate the entire retrieval before writing either file. An altered source
    # requires an explicit new manifest/review; it is never silently substituted.
    for name, data in content.items():
        row = expected[name]
        if len(data) != int(row["bytes"]) or hashlib.md5(data).hexdigest() != row["md5"]:
            raise ValueError(f"Pinned source mismatch: {name}")
    if args.refresh or not all((ROOT / name).exists() for name in expected):
        for name, data in content.items():
            path = ROOT / name
            temporary = path.with_suffix(path.suffix + ".tmp")
            temporary.write_bytes(data)
            temporary.replace(path)
    print(f"{action}: {len(content)} files; all 2,504 samples, one genuine locus")


if __name__ == "__main__":
    main()
