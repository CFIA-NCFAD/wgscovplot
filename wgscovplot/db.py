from typing import List, Dict

from pydantic import BaseModel

from wgscovplot.features import EChartsFeature
from wgscovplot.tools.mosdepth import MosdepthDepthInfo
from wgscovplot.tools.mosdepth.flu import FluMosdepthDepthInfo


class DB(BaseModel):
    samples: List[str]
    low_coverage_threshold: int = 10


class SegmentedGenomeDB(DB):
    """Singleton wgscovplot template DB class containing all data to populate HTML plots

    The single instance of this class will be serialized into a JS object for populating plots in the final HTML output.
    """
    segments: List[str]
    segments_ref_id: Dict[str, Dict[str, str]]
    segments_ref_seq: Dict[str, Dict]
    depths: Dict[str, Dict[str, str]]
    mosdepth_info: Dict[str, Dict[str, FluMosdepthDepthInfo]]
    primer_matches:  Dict[str, Dict[str, List[Dict]]]
    variants: Dict[str, Dict[str, List[Dict]]]


class NonSegmentedGenomeDB(DB):
    """Singleton wgscovplot template DB class containing all data to populate HTML plots

    The single instance of this class will be serialized into a JS object for populating plots in the final HTML output.
    """
    ref_seq: str
    depths: Dict[str, str]
    amplicon_depths: Dict[str, List[Dict]]
    mosdepth_info: Dict[str, MosdepthDepthInfo]
    variants: Dict[str, List[Dict]]
    echart_features: List[EChartsFeature]
