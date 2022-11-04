from collections import defaultdict
from pathlib import Path
from typing import Mapping, List, Dict

import pandas as pd
from Bio import SeqIO

from wgscovplot.flu import read_top_references_table
from wgscovplot.tools.mosdepth import TOP_REFERENCE_PATTERNS, SAMPLE_NAME_CLEANUP, logger, read_mosdepth_bed, \
    depth_array, count_positions, get_interval_coords_bed, get_genome_coverage, get_genome_length, PER_BASE_PATTERNS
from wgscovplot.util import find_file_for_each_sample


def get_segments(sample_top_references: Mapping[str, Path]) -> List:
    segments = []
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        segments = set(segments) | set(df['segment'])
    return sorted(list(segments))


def get_segments_depth(
        basedir: Path,
        sample_top_references: Mapping[str, Path]
) -> Dict[str, Dict[str, List[int]]]:
    out = defaultdict(dict)
    for sample, path in sample_top_references.items():
        df = read_top_references_table(path)
        for row in df.itertuples():
            segment = row.segment
            ncbi_id = row.ncbi_id
            bed_files = list(basedir.glob(f'**/mosdepth/**/{sample}.Segment_{segment}.{ncbi_id}.per-base.bed.gz'))
            if not bed_files:
                logger.error(f'No Mosdepth BED file found for sample "{sample}" and segment "{segment}" ({ncbi_id}).')
                continue
            mosdepth_bed = bed_files[0]
            df_mosdepth = read_mosdepth_bed(mosdepth_bed)
            arr = depth_array(df_mosdepth)
            out[sample][segment] = arr.tolist()
    return dict(out)


def get_segments_ref_seq(
        basedir: Path,
        sample_top_references: Mapping[str, Path]
) -> Dict[str, Dict[str, str]]:
    out = defaultdict(dict)
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        for row in df.itertuples():
            ref_files = basedir.glob(f'**/reference_sequences/**/'
                                     f'{row.sample}.Segment_{row.segment}.{row.ncbi_id}.*')
            for p in ref_files:
                for record in SeqIO.parse(open(p), 'fasta'):
                    out[sample][row.segment] = str(record.seq)
    return dict(out)


def get_segments_ref_id(
        sample_top_references: Mapping[str, Path]
) -> Dict[str, Dict[str, str]]:
    out = defaultdict(dict)
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        for row in df.itertuples():
            out[sample][row.segment] = row.ncbi_id
    return out


def get_flu_info(
        basedir: Path,
        samples: List,
        segments: List,
        sample_top_references: Mapping[str, Path],
        low_coverage_threshold: int = 5
) -> str:
    # TODO: use MosdepthInfo data that should be already be gathered
    headers = [
        'Sample',
        'Segment',
        '# 0 Coverage Positions',
        '0 Coverage Regions',
        f'# Low Coverage Positions (< {low_coverage_threshold}X)',
        f'Low Coverage Regions (< {low_coverage_threshold}X)',
        f'% Genome Coverage >= {low_coverage_threshold}X',
        'Mean Coverage Depth (X)',
        'Median Coverage Depth (X)',
        'Ref Sequence Length (bp)'
    ]
    df_flu_info = pd.DataFrame(
        columns=headers, index=list(range(len(samples) * len(segments)))
    )

    for i, sample in enumerate(samples):
        for j, segment in enumerate(segments):
            df_flu_info.loc[i + i * (len(segments) - 1) + j, ['Sample', 'Segment']] = [sample, segment]
    df_flu_info = df_flu_info.sort_values(
        by=["Sample", "Segment"], ascending=[True, True])
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


def get_low_coverage_regions(
        basedir: Path,
        segments_name: List,
        sample_top_references: Mapping[str, Path],
        low_coverage_threshold: int = 5
) -> Dict[str, Dict[str, str]]:
    out = defaultdict(dict)
    for sample, top_references_path in sample_top_references.items():
        df = read_top_references_table(top_references_path)
        for row in df.itertuples():
            segment = row.segment
            bed_files = basedir.glob(f'**/mosdepth/**/'
                                     f'{sample}.Segment_{segment}.{row.ncbi_id}.per-base.bed.gz')
            for p in bed_files:
                df_mosdepth = read_mosdepth_bed(p)
                out[sample][segment] = get_interval_coords_bed(df_mosdepth, low_coverage_threshold)
    return dict(out)


def get_samples_name(basedir: Path, segment_virus: bool) -> List:
    glob_patterns = TOP_REFERENCE_PATTERNS if segment_virus else PER_BASE_PATTERNS
    sample_beds = find_file_for_each_sample(basedir,
                                            glob_patterns=glob_patterns,
                                            sample_name_cleanup=SAMPLE_NAME_CLEANUP)
    out = []
    for sample, bed_path in sample_beds.items():
        out.append(sample)
    return sorted(out)
