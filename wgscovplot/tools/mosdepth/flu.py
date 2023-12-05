from pathlib import Path
import logging
from typing import List, Tuple, Dict
import base64

import pandas as pd
from Bio import SeqIO

from wgscovplot.tools.mosdepth import SAMPLE_NAME_CLEANUP
from wgscovplot.util import find_file_for_each_sample
from wgscovplot.tools.mosdepth import logger, read_mosdepth_bed, depth_array, count_positions, \
    get_interval_coords_bed, get_genome_coverage, get_genome_length, MosdepthDepthInfo

logger = logging.getLogger(__name__)

TOP_REFERENCE_PATTERNS = [
    '**/reference_sequences/**/*.topsegments.csv',
]


class FluMosdepthDepthInfo(MosdepthDepthInfo):
    segment: str


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


def get_sample_top_references(basedir: Path) -> Dict[str, pd.DataFrame]:
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
    out = {sample: read_top_references_table(top_references_path)
           for sample, top_references_path in sample_top_references.items()}
    return out


def get_segments(sample_top_references: Dict[str, pd.DataFrame]) -> List:
    """
    Return the list of unique segments across all samples
    ['1_PB2', '2_PB1', '3_PA', '4_HA', '5_NP', '6_NA', '7_M', '8_NS']
    """
    segments = []
    for sample, df_top_references in sample_top_references.items():
        segments = set(segments) | set(df_top_references['segment'])
    return sorted(list(segments))


def get_segments_ref_id(
        segments: List,
        sample_top_references: Dict[str, pd.DataFrame]
) -> Dict[str, Dict[str, str]]:
    """
    Return a dictionary of {sample name:{ segment: ref_id}}
    For example:
    'WIN_CFIA_AIV_SAMPLE_96': {'1_PB2': 'MZ171357.1', '2_PB1': 'OQ733010.1', '3_PA': 'OP377327.1', '4_HA': 'OQ733004.1',
                              '5_NP': 'OQ694923.1', '6_NA': 'OQ733005.1', '7_M': 'OQ733006.1', '8_NS': 'OQ694932.1'}
    """
    out = {}
    for sample, df_top_references in sample_top_references.items():
        out[sample] = {}
        for segment in segments:
            # some samples may not have enough segments, if sample, segment has no info so set empty
            out[sample][segment] = ""
        for row in df_top_references.itertuples():
            out[sample][row.segment] = row.ref_id
    return out


def get_segments_ref_seq(
        basedir: Path,
        segments: List,
        sample_top_references: Dict[str, pd.DataFrame]
) -> Dict[str, Dict[str, str]]:
    """
    Return a dictionary of {sample name:{ segment: ref_seq}}
    For example:
    'WIN_CFIA_AIV_SAMPLE_96': {'1_PB2': 'ACGCGCAGCAGCG', '2_PB1': 'ACGCGCAGCAGCGGC', '3_PA': 'ACGCGCAGCAGCG', '4_HA': 'ACGCGCAGCAGCG',
                              '5_NP': 'ACGCGCAGCAGCG', '6_NA': 'ACGCGCAGCAGCG', '7_M': 'ACGCGCAGCAGCG', '8_NS': 'ACGCGCAGCAGCG'}
    """
    out = {}
    for sample, df_top_references in sample_top_references.items():
        out[sample] = {}
        for segment in segments:
            out[sample][segment] = ""  # if sample, segment has no info so set empty
        for row in df_top_references.itertuples():
            segment = row.segment
            ref_id = row.ref_id
            ref_files = list(basedir.glob(f'**/reference_sequences/**/'
                                          f'{row.sample}.Segment_{segment}.{ref_id}.*'))
            if not ref_files:
                logger.error(f'No Reference file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
                continue
            ref_file = ref_files[0]
            for record in SeqIO.parse(open(ref_file), 'fasta'):
                out[sample][segment] = str(record.seq)
                #out[sample][segment] = zlib.compress(str(record.seq).encode())
                # Test size after and before compress:
                # print(sys.getsizeof(zlib.compress(str(record.seq).encode())), sys.getsizeof(str(record.seq)))
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
    out = {}
    for sample, df_top_references in sample_top_references.items():
        out[sample] = {}
        for segment in segments:
            out[sample][segment] = ""  # if sample, segment has no info so set empty
        for row in df_top_references.itertuples():
            segment = row.segment
            ref_id = row.ref_id
            consensus_files = list(basedir.glob(f'**/consensus/bcftools/**/'
                                                f'{row.sample}.Segment_{segment}.{ref_id}.*'))
            if not consensus_files:
                logger.error(f'No Consensus file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
                continue
            consensus_file = consensus_files[0]
            for record in SeqIO.parse(open(consensus_file), 'fasta'):
                out[sample][segment] = str(record.seq)
    return out


def get_flu_mosdepth_info(
        basedir: Path,
        segments: List,
        sample_top_references: Dict[str, pd.DataFrame],
        low_coverage_threshold: int = 5
) -> Tuple[Dict[str, Dict[str, MosdepthDepthInfo]], Dict[str, Dict[str, str]]]:
    out = {}
    sample_depths = {}
    for sample, df_top_references in sample_top_references.items():
        sample_depths[sample] = {}
        out[sample] = {}
        for segment in segments:
            sample_depths[sample][segment] = ""  # if sample, segment has no info so set empty
            out[sample][segment] = FluMosdepthDepthInfo(sample=sample,
                                                        segment=segment,
                                                        low_coverage_threshold=low_coverage_threshold)
        for row in df_top_references.itertuples():
            segment = row.segment
            ref_id = row.ref_id
            bed_files = list(basedir.glob(f'**/mosdepth/**/{sample}.Segment_{segment}.{ref_id}.per-base.bed.gz'))
            if not bed_files:
                logger.error(f'No Mosdepth BED file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
                continue
            mosdepth_bed = bed_files[0]
            df_mosdepth = read_mosdepth_bed(mosdepth_bed)
            arr = depth_array(df_mosdepth)
            mean_cov = arr.mean()
            median_cov = pd.Series(arr).median()
            depth_info = FluMosdepthDepthInfo(sample=sample,
                                              segment=segment,
                                              low_coverage_threshold=low_coverage_threshold,
                                              n_low_coverage=count_positions(df_mosdepth[df_mosdepth.depth <
                                                                                         low_coverage_threshold]),
                                              n_zero_coverage=count_positions(df_mosdepth[df_mosdepth.depth == 0]),
                                              zero_coverage_coords=get_interval_coords_bed(df_mosdepth),
                                              low_coverage_coords=get_interval_coords_bed(df_mosdepth,
                                                                                          low_coverage_threshold),
                                              genome_coverage=get_genome_coverage(df_mosdepth,
                                                                                  low_coverage_threshold),
                                              mean_coverage=mean_cov,
                                              median_coverage=median_cov,
                                              ref_seq_length=get_genome_length(df_mosdepth),
                                              max_depth=arr.max())
            out[sample][segment] = depth_info
            arr[arr == 0] = 1E-7  # assign values for zero depth position so that it can be plotted in log scale mode
            sample_depths[sample][segment] = base64.b64encode(arr).decode('utf-8')
    return out, sample_depths
