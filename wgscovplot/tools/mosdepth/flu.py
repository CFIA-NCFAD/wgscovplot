import logging
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Tuple

import numpy as np
import pandas as pd

from wgscovplot.flu import SampleSegmentRef
from wgscovplot.tools.mosdepth import (
    MosdepthDepthInfo,
    count_positions,
    depth_array,
    get_genome_coverage,
    get_genome_length,
    get_interval_coords_bed,
    read_mosdepth_bed,
)

logger = logging.getLogger(__name__)


class FluMosdepthDepthInfo(MosdepthDepthInfo):
    segment: str


def get_flu_mosdepth_info(
    basedir: Path,
    sample_top_references: List[SampleSegmentRef],
    low_coverage_threshold: int = 5,
) -> Tuple[Dict[str, Dict[str, FluMosdepthDepthInfo]], Dict[str, Dict[str, np.ndarray]]]:
    sample_segment_depths: DefaultDict[str, Dict] = defaultdict(dict)
    sample_segment_depth_info: DefaultDict[str, Dict] = defaultdict(dict)
    for item in sample_top_references:
        segment = item.segment
        ref_id = item.ref_id
        sample = item.sample
        bed_files = list(basedir.glob(f"**/mosdepth/**/{sample}.Segment_{segment}.{ref_id}.per-base.bed.gz"))
        if not bed_files:
            logger.error(f'No Mosdepth BED file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
            continue
        mosdepth_bed = bed_files[0]
        df_mosdepth = read_mosdepth_bed(mosdepth_bed)
        arr = depth_array(df_mosdepth)
        mean_cov = arr.mean()
        median_cov = pd.Series(arr).median()
        depth_info = FluMosdepthDepthInfo(
            sample=sample,
            segment=segment,
            low_coverage_threshold=low_coverage_threshold,
            n_low_coverage=count_positions(df_mosdepth[df_mosdepth.depth < low_coverage_threshold]),
            n_zero_coverage=count_positions(df_mosdepth[df_mosdepth.depth == 0]),
            zero_coverage_coords=get_interval_coords_bed(df_mosdepth),
            low_coverage_coords=get_interval_coords_bed(df_mosdepth, low_coverage_threshold),
            genome_coverage=get_genome_coverage(df_mosdepth, low_coverage_threshold),
            mean_coverage=mean_cov,
            median_coverage=median_cov,
            ref_seq_length=get_genome_length(df_mosdepth),
            max_depth=arr.max(),
        )
        sample_segment_depth_info[sample][segment] = depth_info
        sample_segment_depths[sample][segment] = arr
    return dict(sample_segment_depth_info), dict(sample_segment_depths)
