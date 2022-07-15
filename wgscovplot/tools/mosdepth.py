import logging
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO

import numpy as np
import pandas as pd
from pydantic import BaseModel

from wgscovplot.colors import AmpliconColour
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
    if '.tsv' in Path(p).suffixes:
        df_temp = pd.read_table(p, header=None, names=['sample_name', 'reference', 'pos', 'depth'])
        coverted_df = pd.DataFrame(columns=['genome', 'start_idx', 'end_idx', 'depth'])
        coverted_df['genome'] = df_temp['reference']
        coverted_df['start_idx'] = df_temp['pos'] - 1
        coverted_df['end_idx'] = df_temp['pos']
        coverted_df['depth'] = df_temp['depth']
        return coverted_df
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


def get_refseq_name(basedir: Path) -> str:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    refseq_name = ''
    for sample, bed_path in sample_beds.items():
        df = read_mosdepth_bed(bed_path)
        refseq_name = df['genome'][0]
        break
    return refseq_name


def get_depth(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=PER_BASE_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    for sample, bed_path in sample_beds.items():
        df = read_mosdepth_bed(bed_path)
        arr = depth_array(df)
        arr[arr == 0] = 1E-20
        out[sample] = arr.tolist()
    return out


def get_samples_name(basedir: Path, segment_virus: bool) -> List:
    glob_patterns = TOP_REFERENCE_PATTERNS if segment_virus else PER_BASE_PATTERNS
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=glob_patterns,
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
    try:
        for sample, bed_path in sample_beds.items():
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
    except:
        logging.warning('No Region Amplicon Depth Found')
        return {}


def get_region_amplicon(basedir: Path) -> Dict[str, List]:
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=REGIONS_PATTERNS,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    try:
        for sample, bed_path in sample_beds.items():
            df_amplicon = pd.read_table(bed_path,
                                        names=['reference', 'start', 'end', 'amplicon', 'depth'],
                                        header=None)
            out = {row.amplicon: [row.start, row.end] for row in df_amplicon.itertuples()}
            break
        return out
    except:
        logging.warning('No Region Amplicon Found')
        return {}


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


def read_top_references_table(basedir: Path) -> pd.DataFrame:
    return pd.read_csv(basedir, sep=',', header=0, names=['sample', 'segment', 'ncbi_id',
                                                          'blastn_bitscore', 'ref_sequence_id'])


def get_segments_depth(basedir: Path) -> Dict[str, Dict[str, List]]:
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    segments_name = get_segments_name(basedir)
    for sample, top_references_path in sample_top_references.items():
        out[sample] = {}
        df = read_top_references_table(top_references_path)
        for segment in segments_name:
            out[sample][segment] = []
        for row in df.itertuples():
            bed_files = basedir.glob(f'**/mosdepth/**/'
                                     f'{row.sample}.Segment_{row.segment}.{row.ncbi_id}.per-base.bed.gz')
            for p in bed_files:
                df_mosdepth = read_mosdepth_bed(p)
                arr = depth_array(df_mosdepth)
                arr[arr == 0] = 1E-20
                out[sample][row.segment] = arr.tolist()
    return out


def get_segments_references(basedir: Path) -> Dict[str, List]:
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    segments_name = get_segments_name(basedir)
    for sample, top_references_path in sample_top_references.items():
        out[sample] = {}
        df = read_top_references_table(top_references_path)
        for segment in segments_name:
            out[sample][segment] = ''
        for row in df.itertuples():
            ref_files = basedir.glob(f'**/reference_sequences/**/'
                                     f'{row.sample}.Segment_{row.segment}.{row.ncbi_id}.*')
            for p in ref_files:
                for record in SeqIO.parse(open(p), 'fasta'):
                    out[sample][row.segment] = str(record.seq)
    return out


def get_segments_ref_id(basedir: Path) -> Dict[str, Dict[str, str]]:
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    segments_name = get_segments_name(basedir)
    for sample, top_references_path in sample_top_references.items():
        out[sample] = {}
        df = read_top_references_table(top_references_path)
        for segment in segments_name:
            out[sample][segment] = ''
        for row in df.itertuples():
            out[sample][row.segment] = row.ncbi_id
    return out


def get_segments_name(basedir: Path) -> List:
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    segments = []
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        segments = set(segments) | set(df['segment'])
    return sorted(list(segments))


def get_flu_info(basedir: Path, samples: List, segments: List, low_coverage_threshold: int = 5) -> str:
    headers = ['Sample', 'Segment', '# 0 Coverage Positions', '0 Coverage Regions',
               '# Low Coverage Positions (< ' + str(low_coverage_threshold) + 'X)',
               'Low Coverage Regions (< ' + str(low_coverage_threshold) + 'X)',
               '% Genome Coverage >= ' + str(low_coverage_threshold) + 'X',
               'Mean Coverage Depth (X)',
               'Median Coverage Depth (X)', 'Ref Sequence Length (bp)']
    df_flu_info = pd.DataFrame(columns=headers, index=list(range(0, len(samples) * len(segments))))
    for i, sample in enumerate(samples):
        for j, segment in enumerate(segments):
            df_flu_info.loc[i + i * (len(segments) - 1) + j, ['Sample', 'Segment']] = [sample, segment]
    df_flu_info = df_flu_info.sort_values(
        by=["Sample", "Segment"], ascending=[True, True])
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        for row in df.itertuples():
            bed_files = basedir.glob(f'**/mosdepth/**/'
                                     f'{row.sample}.Segment_{row.segment}.{row.ncbi_id}.per-base.bed.gz')
            for p in bed_files:
                df_mosdepth = read_mosdepth_bed(p)
                arr = depth_array(df_mosdepth)
                mean_cov = arr.mean()
                median_cov = pd.Series(arr).median()
                depth_info = [row.sample,
                              row.segment,
                              count_positions(df_mosdepth[df_mosdepth.depth == 0]),
                              get_interval_coords_bed(df_mosdepth),
                              count_positions(df_mosdepth[df_mosdepth.depth < low_coverage_threshold]),
                              get_interval_coords_bed(df_mosdepth, low_coverage_threshold),
                              get_genome_coverage(df_mosdepth, low_coverage_threshold) * 100,
                              mean_cov,
                              median_cov,
                              get_genome_length(df_mosdepth)]
                df_flu_info.loc[
                    (df_flu_info['Sample'] == row.sample) & (df_flu_info['Segment'] == row.segment)] = depth_info
    df_flu_info.fillna('No result reported', inplace=True)
    df_flu_info.index = df_flu_info.index + 1
    return df_flu_info.to_html(classes="table table-striped table-hover table-bordered table-responsive-md",
                               float_format=lambda x: f'{x:0.2f}', justify="left", table_id="summary-coverage-stat")


def get_low_coverage_regions(basedir: Path, low_coverage_threshold: int = 5) -> Dict[str, Dict[str, str]]:
    sample_top_references = find_file_for_each_sample(basedir,
                                                      glob_patterns=TOP_REFERENCE_PATTERNS,
                                                      sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = {}
    segments_name = get_segments_name(basedir)
    for sample, top_references_path in sample_top_references.items():
        out[sample] = {}
        df = read_top_references_table(top_references_path)
        for segment in segments_name:
            out[sample][segment] = ''
        for row in df.itertuples():
            bed_files = basedir.glob(f'**/mosdepth/**/'
                                     f'{row.sample}.Segment_{row.segment}.{row.ncbi_id}.per-base.bed.gz')
            for p in bed_files:
                df_mosdepth = read_mosdepth_bed(p)
                out[sample][row.segment] = get_interval_coords_bed(df_mosdepth, low_coverage_threshold)
    return out
