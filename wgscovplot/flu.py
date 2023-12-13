from collections import defaultdict
from pathlib import Path
import pandas as pd
import logging
from Bio.SeqIO.FastaIO import SimpleFastaParser
from typing import Dict, List

from wgscovplot.util import find_file_for_each_sample
from wgscovplot.tools.mosdepth.flu import TOP_REFERENCE_PATTERNS, SampleSegmentRef
from wgscovplot.tools.mosdepth import SAMPLE_NAME_CLEANUP

logger = logging.getLogger(__name__)


def get_sample_top_references(basedir: Path) -> List[SampleSegmentRef]:
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
    df = pd.concat([read_top_references_table(p) for _, p in sample_top_references.items()])
    return [
        SampleSegmentRef(
            sample=r.sample,
            segment=r.segment,
            ref_id=r.ref_id,
        ) for r in df.itertuples()]


def get_segments(sample_top_references: List[SampleSegmentRef]) -> List:
    """
    Return the list of unique segments across all samples
    ['1_PB2', '2_PB1', '3_PA', '4_HA', '5_NP', '6_NA', '7_M', '8_NS']
    """
    segments = []
    for item in sample_top_references:
        if item.segment not in segments:
            segments.append(item.segment)
    return sorted(list(segments))


def read_top_references_table(basedir: Path) -> pd.DataFrame:
    """
    Parse sample.topsegments.csv into data frame with five columns:
            'sample',
            'segment',
            'ref_id',
            'blastn_bitscore',
            'ref_sequence_id',
    """
    return pd.read_csv(
        basedir,
        header=0,
        names=[
            'sample',
            'segment',
            'ref_id',
            'blastn_bitscore',
            'ref_sequence_id',
        ]
    )


def get_segments_ref_id(
        segments: List,
        sample_top_references: List[SampleSegmentRef]
) -> Dict[str, Dict[str, str]]:
    """
    Return a dictionary of {sample name:{ segment: ref_id}}
    For example:
    'WIN_CFIA_AIV_SAMPLE_96': {'1_PB2': 'MZ171357.1', '2_PB1': 'OQ733010.1', '3_PA': 'OP377327.1', '4_HA': 'OQ733004.1',
                              '5_NP': 'OQ694923.1', '6_NA': 'OQ733005.1', '7_M': 'OQ733006.1', '8_NS': 'OQ694932.1'}
    """

    out = defaultdict(dict)
    for items in sample_top_references:
        out[items.sample][items.segment] = items.ref_id
    return out


def get_segments_ref_seq(
        basedir: Path,
        segments: List,
        sample_top_references: SampleSegmentRef
) -> Dict[str, Dict[str, str]]:
    """
    Return a dictionary of {sample name:{ segment: ref_seq}}
    For example:
    'WIN_CFIA_AIV_SAMPLE_96': {'1_PB2': 'ACGCGCAGCAGCG', '2_PB1': 'ACGCGCAGCAGCGGC', '3_PA': 'ACGCGCAGCAGCG', '4_HA': 'ACGCGCAGCAGCG',
                              '5_NP': 'ACGCGCAGCAGCG', '6_NA': 'ACGCGCAGCAGCG', '7_M': 'ACGCGCAGCAGCG', '8_NS': 'ACGCGCAGCAGCG'}
    """
    out = defaultdict(dict)
    for items in sample_top_references:
        segment = items.segment
        ref_id = items.ref_id
        sample = items.sample
        ref_files = list(basedir.glob(f'**/reference_sequences/**/'
                                      f'{sample}.Segment_{segment}.{ref_id}.*'))
        if not ref_files:
            logger.error(f'No Reference file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
            continue
        ref_file = ref_files[0]
        with open(ref_file) as handle:
            for _, seq in SimpleFastaParser(handle):
                out[sample][segment] = str(seq)

    return out


def get_segments_consensus_seq(
        basedir: Path,
        segments: List,
        sample_top_references: Dict[str, pd.DataFrame]
) -> Dict[str, Dict[str, str]]:
    """
    Return a dictionary of {sample name:{ segment: consensus_sequence}}
    For example:
    'WIN_CFIA_AIV_SAMPLE_96': {'1_PB2': 'ACGCGCAGCAGCG', '2_PB1': 'ACGCGCAGCAGCGGC', '3_PA': 'ACGCGCAGCAGCG', '4_HA': 'ACGCGCAGCAGCG',
                              '5_NP': 'ACGCGCAGCAGCG', '6_NA': 'ACGCGCAGCAGCG', '7_M': 'ACGCGCAGCAGCG', '8_NS': 'ACGCGCAGCAGCG'}
    """
    out = defaultdict(dict)
    for items in sample_top_references:
        segment = items.segment
        ref_id = items.ref_id
        sample = items.sample
        consensus_files = list(basedir.glob(f'**/consensus/bcftools/**/'
                                            f'{sample}.Segment_{segment}.{ref_id}.*'))
        if not consensus_files:
            logger.error(f'No Consensus file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
            continue
        consensus_file = consensus_files[0]
        with open(consensus_file) as handle:
            for _, seq in SimpleFastaParser(handle):
                out[sample][segment] = str(seq)
    return out
