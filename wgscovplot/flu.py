import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pydantic import BaseModel

from wgscovplot.tools.mosdepth import SAMPLE_NAME_CLEANUP
from wgscovplot.util import find_file_for_each_sample, select_most_recent_file

logger = logging.getLogger(__name__)

SEGMENTS = ["1_PB2", "2_PB1", "3_PA", "4_HA", "5_NP", "6_NA", "7_M", "8_NS"]


class SampleSegmentRef(BaseModel):
    sample: str
    segment: str
    ref_id: str


def get_sample_top_references(basedir: Path) -> list[SampleSegmentRef]:
    """
    Find files sample.topsegments.csv for each sample
    File is located in reference_sequences/sample/sample.topsegments.csv
    Return a dictionary { sample_name: data_frame_with_five_columns }
    """
    sample_top_references = find_file_for_each_sample(
        basedir,
        glob_patterns=TOP_REFERENCE_PATTERNS,
        sample_name_cleanup=SAMPLE_NAME_CLEANUP,
    )
    if not sample_top_references:
        return []
    df = pd.concat([read_nf_flu_topsegments_csv(p) for _, p in sample_top_references.items()])
    return [
        SampleSegmentRef(
            sample=r.sample,
            segment=r.segment,
            ref_id=r.ref_id,
        )
        for r in df.itertuples()
    ]


def read_nf_flu_topsegments_csv(top_segments_csv: Path) -> pd.DataFrame:
    return pd.read_csv(
        top_segments_csv,
        header=0,
        names=[
            "sample",
            "segment",
            "ref_id",
            "blastn_bitscore",
            "ref_sequence_id",
        ],
    )


def get_segments_ref_id(sample_top_references: list[SampleSegmentRef]) -> dict[str, dict[str, str]]:
    """Get reference id for each segment of each sample"""
    out: defaultdict[str, dict[str, str]] = defaultdict(dict)
    for items in sample_top_references:
        out[items.sample][items.segment] = items.ref_id
    return dict(out)


def get_sample_segment_seqs(
    basedir: Path,
    sample_top_references: list[SampleSegmentRef],
) -> dict[str, dict[str, str]]:
    """Get sequence for each segment of each sample

    This function assumes that the nf-flu pipeline was used, so the consensus sequence would have the filename format:
    {sample}.Segment_{segment}.{ref_id}.bcftools.consensus.fasta
    """
    out: defaultdict[str, dict[str, str]] = defaultdict(dict)
    # Find the most recently modified FASTA file for each sample for each segment
    # Only do the glob search once rather than once for every sample and segment
    filename_to_fasta_file = find_file_for_each_sample(
        basedir,
        glob_patterns=["**/*.fasta"],
        sample_name_cleanup=[".reference", ".fasta"],
        single_entry_selector_func=select_most_recent_file,
    )
    logger.info(f"Found {len(filename_to_fasta_file)} FASTA files in '{basedir}'")
    for ssr in sample_top_references:
        segment = ssr.segment
        ref_id = ssr.ref_id
        sample = ssr.sample
        filename = f"{sample}.Segment_{segment}.{ref_id}"
        fasta_file = filename_to_fasta_file.get(filename)
        if not fasta_file:
            logger.error(f'No FASTA file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
            continue
        with open(fasta_file) as handle:
            for _, seq in SimpleFastaParser(handle):
                out[sample][segment] = seq
    if not out:
        msg = "No FASTA files found for any sample and segment"
        raise FileNotFoundError(msg)
    return dict(out)


TOP_REFERENCE_PATTERNS = [
    "**/*.topsegments.csv",
]
