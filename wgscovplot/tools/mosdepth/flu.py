import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from wgscovplot.flu import SampleSegmentRef
from wgscovplot.tools.mosdepth import (
    MosdepthDepthInfo,
    depth_array,
    get_genome_coverage,
    get_interval_coords,
    read_mosdepth_bed,
)

logger = logging.getLogger(__name__)


class FluMosdepthDepthInfo(MosdepthDepthInfo):
    segment: str


def get_flu_mosdepth_info(
    basedir: Path,
    sample_top_references: list[SampleSegmentRef],
    low_coverage_threshold: int = 5,
) -> tuple[dict[str, dict[str, FluMosdepthDepthInfo]], dict[str, dict[str, np.ndarray]]]:
    sample_segment_depths: defaultdict[str, dict] = defaultdict(dict)
    sample_segment_depth_info: defaultdict[str, dict] = defaultdict(dict)
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
            n_low_coverage=(arr < low_coverage_threshold).sum(),
            n_zero_coverage=(arr == 0).sum(),
            zero_coverage_coords=get_interval_coords(arr),
            low_coverage_coords=get_interval_coords(arr, low_coverage_threshold),
            genome_coverage=get_genome_coverage(arr, low_coverage_threshold),
            mean_coverage=mean_cov,
            median_coverage=median_cov,
            ref_seq_length=len(arr),
            max_depth=arr.max(initial=0),
        )
        sample_segment_depth_info[sample][segment] = depth_info
        sample_segment_depths[sample][segment] = arr
    return dict(sample_segment_depth_info), dict(sample_segment_depths)
