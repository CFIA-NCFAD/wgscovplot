import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from pydantic import BaseModel

from wgscovplot.colors import AmpliconColour
from wgscovplot.util import find_file_for_each_sample

logger = logging.getLogger(__name__)

SAMPLE_NAME_CLEANUP = [
    '.genome.per-base.bed.gz',
    '.per-base.bed.gz',
    '.genome.regions.bed.gz',
    '.regions.bed.gz',
    '.trim',
    '.mkD',

]

PER_BASE_PATTERNS = [
    '**/mosdepth/**/*.genome.per-base.bed.gz',
    '**/mosdepth/**/*.per-base.bed.gz',
]

REGIONS_PATTERNS = [
    '**/mosdepth/**/*.genome.regions.bed.gz',
    '**/mosdepth/**/*.regions.bed.gz'
]


class MosdepthDepthInfo(BaseModel):
    sample: str
    n_zero_coverage: int
    zero_coverage_coords: str
    low_coverage_threshold: int = 5
    n_low_coverage: int
    low_coverage_coords: str
    genome_coverage: float
    mean_coverage: float
    median_coverage: int
    ref_seq_length: int


def read_mosdepth_bed(p: Path) -> pd.DataFrame:
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
    arr = np.zeros(df.end_idx.max())
    for row in df.itertuples():
        arr[row.start_idx:row.end_idx] = row.depth
    return arr


def read_regular_depths(fpath: Path) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['sample_name', 'reference', 'pos', 'depth'],
                       header=None)
    return df


def get_depth(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    for sample, bed_path in sample_beds.items():
        logger.info(f'Sample "{sample}" has depth file "{bed_path}"')
        df = read_mosdepth_bed(bed_path)
        arr = depth_array(df)
        arr[arr == 0] = 1E-20
        out[sample] = arr.tolist()
    return out


def get_samples_name(basedir: Path) -> List:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = []
    for sample, bed_path in sample_beds.items():
        out.append(sample)
    return sorted(out)


def get_depth_amplicon(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=REGIONS_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    for sample, bed_path in sample_beds.items():
        logger.info(f'Sample "{sample}" has region depth file "{bed_path}"')
        df = read_mosdepth_region_bed(bed_path)
        out[sample] = []
        for row in df.itertuples():
            pool_id = int(row.amplicon.split('_')[-1])
            if pool_id % 2:  # pool 2
                out[sample].append(dict(
                    value=[row.start_idx, row.end_idx, row.depth, row.amplicon],
                    itemStyle={"color": AmpliconColour.pool2.value}))
            else:  # pool 1
                out[sample].append(dict(
                    value=[row.start_idx, row.end_idx, row.depth, row.amplicon],
                    itemStyle={"color": AmpliconColour.pool1.value}))
    return out


def get_region_amplicon(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=REGIONS_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    for sample, bed_path in sample_beds.items():
        df_amplicon = pd.read_table(bed_path,
                                    names=['reference', 'start', 'end', 'amplicon', 'depth'],
                                    header=None)
        out = {row.amplicon: [row.start, row.end] for row in df_amplicon.itertuples()}
        break
    return out


def get_info(basedir: Path, low_coverage_threshold: int = 5) -> Dict[str, MosdepthDepthInfo]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    for sample, bed_path in sample_beds.items():
        df = read_mosdepth_bed(bed_path)
        arr = depth_array(df)
        mean_cov = arr.mean()
        median_cov = pd.Series(arr).median()
        depth_info = MosdepthDepthInfo(sample=sample,
                                       low_coverage_threshold=low_coverage_threshold,
                                       n_low_coverage=count_positions(df[df.depth < low_coverage_threshold]),
                                       n_zero_coverage=count_positions(df[df.depth == 0]),
                                       zero_coverage_coords=get_interval_coords_bed(df),
                                       low_coverage_coords=get_interval_coords_bed(df, low_coverage_threshold),
                                       genome_coverage=get_genome_coverage(df, low_coverage_threshold),
                                       mean_coverage=mean_cov,
                                       median_coverage=median_cov,
                                       ref_seq_length=get_genome_length(df))
        out[sample] = depth_info
    return out
