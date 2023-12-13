from pathlib import Path
from collections import defaultdict
import logging
from typing import List, Tuple, Dict
import base64
import pandas as pd
from pydantic.main import BaseModel

from wgscovplot.tools.mosdepth import read_mosdepth_bed, depth_array, count_positions, \
    get_interval_coords_bed, get_genome_coverage, get_genome_length, MosdepthDepthInfo


logger = logging.getLogger(__name__)

TOP_REFERENCE_PATTERNS = [
    '**/*.topsegments.csv',
]


class FluMosdepthDepthInfo(MosdepthDepthInfo):
    segment: str


class SampleSegmentRef(BaseModel):
    sample: str
    segment: str
    ref_id: str


def get_flu_mosdepth_info(
        basedir: Path,
        segments: List,
        sample_top_references: SampleSegmentRef,
        low_coverage_threshold: int = 5
) -> Tuple[Dict[str, Dict[str, MosdepthDepthInfo]], Dict[str, Dict[str, str]]]:
    sample_depths = defaultdict(dict)
    mosdepth_info = defaultdict(dict)
    for items in sample_top_references:
        segment = items.segment
        ref_id = items.ref_id
        sample = items.sample
        bed_files = list(basedir.glob(f'**/mosdepth/**/{sample}.Segment_{segment}.{ref_id}.per-base.bed.gz'))
        if not bed_files:
            logger.error(f'No Mosdepth BED file found for sample "{sample}" and segment "{segment}" ({ref_id}).')
            continue
        mosdepth_bed = bed_files[0]
        df_mosdepth = read_mosdepth_bed(mosdepth_bed)
        arr = depth_array(df_mosdepth)
        mean_cov = "{:.2f}".format(arr.mean())
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
                                          genome_coverage="{:.2%}".format(get_genome_coverage(df_mosdepth,
                                                                                              low_coverage_threshold)),
                                          mean_coverage=mean_cov,
                                          median_coverage=median_cov,
                                          ref_seq_length=get_genome_length(df_mosdepth),
                                          max_depth=arr.max())
        mosdepth_info[sample][segment] = depth_info
        sample_depths[sample][segment] = base64.b64encode(arr).decode('utf-8')
    return mosdepth_info, sample_depths
