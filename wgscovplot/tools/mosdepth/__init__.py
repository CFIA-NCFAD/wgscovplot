import base64
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from pydantic import BaseModel

from wgscovplot.colors import AmpliconColour
from wgscovplot.features import Feature
from wgscovplot.util import find_file_for_each_sample

logger = logging.getLogger(__name__)

SAMPLE_NAME_CLEANUP = [
    '.genome.per-base.bed.gz',
    '.amplicon.per-base.bed.gz',
    '.per-base.bed.gz',
    '.genome.regions.bed.gz',
    '.amplicon.regions.bed.gz',
    '.regions.bed.gz',
    '-depths.tsv',
    '.trim',
    '.mkD',
    '.topsegments.csv'
]

PER_BASE_PATTERNS = [
    '**/mosdepth/**/*.genome.per-base.bed.gz',
    '**/mosdepth/**/*.per-base.bed.gz',
    '**/mosdepth/**/*-depths.tsv',
]

TOP_REFERENCE_PATTERNS = [
    '**/reference_sequences/**/*.topsegments.csv',
]

REGIONS_PATTERNS = [
    '**/mosdepth/**/*.amplicon.regions.bed.gz',
    '**/mosdepth/**/*.regions.bed.gz'
]


class MosdepthDepthInfo(BaseModel):
    sample: str
    n_zero_coverage: int = 0
    zero_coverage_coords: str = ""
    low_coverage_threshold: int = 10
    n_low_coverage: int = 0
    low_coverage_coords: str = ""
    genome_coverage: float = 0
    mean_coverage: float = 0
    median_coverage: int = 0
    ref_seq_length: int = 0
    max_depth: int = 0


class FluMosdepthDepthInfo(MosdepthDepthInfo):
    segment: str


def read_mosdepth_bed(p: Path) -> pd.DataFrame:
    if '.tsv' in Path(p).suffixes:  # add this code as temporarily support depths.tsv file
        df_temp = pd.read_table(p, header=None, names=['sample_name', 'reference', 'pos', 'depth'])
        converted_df = pd.DataFrame(columns=['genome', 'start_idx', 'end_idx', 'depth'])
        converted_df['genome'] = df_temp['reference']
        converted_df['start_idx'] = df_temp['pos'] - 1
        converted_df['end_idx'] = df_temp['pos']
        converted_df['depth'] = df_temp['depth']
        return converted_df
    return pd.read_table(p, header=None, names=['genome', 'start_idx', 'end_idx', 'depth'])


def read_mosdepth_region_bed(p: Path) -> pd.DataFrame:
    return pd.read_table(p, header=None, names=['genome', 'start_idx', 'end_idx', 'amplicon', 'depth'])


def get_interval_coords_bed(df: pd.DataFrame, threshold: int = 0) -> str:
    mask = df.depth == 0 if threshold == 0 else df.depth < threshold
    df_below = df[mask]
    start_pos, end_pos = df_below.start_idx, df_below.end_idx
    coords = []
    for x, y in zip(start_pos, end_pos):
        if coords:
            last = coords[-1][-1]
            if x == last + 1:
                coords[-1].append(x)
                coords[-1].append(y - 1)
            else:
                coords.append([x, y - 1])
        else:
            coords.append([x, y - 1])
    return '; '.join([f'{xs[0] + 1}-{xs[-1] + 1}' if xs[0] != xs[-1] else f'{xs[0] + 1}' for xs in coords])


def count_positions(df: pd.DataFrame) -> int:
    return (df.end_idx - df.start_idx).sum()


def get_genome_coverage(df: pd.DataFrame, low_coverage_threshold: int = 5) -> float:
    genome_length = get_genome_length(df)
    return 1.0 - (count_positions(df[df.depth < low_coverage_threshold]) / genome_length)


def get_genome_length(df):
    return df.end_idx.max()


def depth_array(df: pd.DataFrame) -> np.ndarray:
    arr = np.zeros(df.end_idx.max(), dtype=np.float32)
    for row in df.itertuples():
        arr[row.start_idx:row.end_idx] = row.depth
    return arr


def get_refseq_id(basedir: Path) -> str:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    refseq_name = ''
    for sample, bed_path in sample_beds.items():
        df = read_mosdepth_bed(bed_path)
        refseq_name = df['genome'][0]
        break
    return refseq_name


def get_samples_name(basedir: Path, segment_virus: bool) -> List:
    glob_patterns = TOP_REFERENCE_PATTERNS if segment_virus else PER_BASE_PATTERNS
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=glob_patterns,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = []
    for sample, bed_path in sample_beds.items():
        out.append(sample)
    return sorted(out)


def get_amplicon_depths(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=REGIONS_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = defaultdict(list)
    sample = None
    try:
        for sample, bed_path in sample_beds.items():
            df = read_mosdepth_region_bed(bed_path)
            for row in df.itertuples():
                pool_id = int(row.amplicon.split('_')[-1])
                color = AmpliconColour.pool2 if pool_id % 2 == 0 else AmpliconColour.pool1
                out[sample].append(
                    dict(
                        value=[
                            row.start_idx,
                            row.end_idx,
                            row.depth,
                            row.amplicon
                        ],
                        itemStyle=dict(color=color),
                    )
                )
        return dict(out)
    except Exception as e:
        logger.error(e, exc_info=True)
        logger.warning(f'{sample} No Region Amplicon Depth Found')
        return {}


def get_region_amplicon(basedir: Path) -> List[Feature]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=REGIONS_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    try:
        for sample, bed_path in sample_beds.items():
            df_amplicon = pd.read_table(bed_path,
                                        names=['reference', 'start', 'end', 'amplicon', 'depth'],
                                        header=None)
            return [
                Feature(
                    name=row.amplicon,
                    start=row.start,
                    end=row.end,
                ) for row in df_amplicon.itertuples()]
    except Warning:
        logger.warning('No Region Amplicon Found')
    return []


def get_info(
        basedir: Path,
        low_coverage_threshold: int = 5
) -> Tuple[Dict[str, MosdepthDepthInfo], Dict[str, np.ndarray]]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    sample_depths = {}
    for sample, bed_path in sample_beds.items():
        df_mosdepth = read_mosdepth_bed(bed_path)
        arr = depth_array(df_mosdepth)
        mean_cov = arr.mean()
        median_cov = pd.Series(arr).median()
        depth_info = MosdepthDepthInfo(sample=sample,
                                       low_coverage_threshold=low_coverage_threshold,
                                       n_low_coverage=count_positions(df_mosdepth[df_mosdepth.depth < low_coverage_threshold]),
                                       n_zero_coverage=count_positions(df_mosdepth[df_mosdepth.depth == 0]),
                                       zero_coverage_coords=get_interval_coords_bed(df_mosdepth),
                                       low_coverage_coords=get_interval_coords_bed(df_mosdepth, low_coverage_threshold),
                                       genome_coverage=get_genome_coverage(df_mosdepth, low_coverage_threshold),
                                       mean_coverage=mean_cov,
                                       median_coverage=median_cov,
                                       ref_seq_length=get_genome_length(df_mosdepth),
                                       max_depth=arr.max())
        out[sample] = depth_info
        arr[arr == 0] = 1E-5
        sample_depths[sample] = base64.b64encode(arr).decode('utf-8')
    return out, sample_depths
