from typing import List, Dict

from pydantic import BaseModel

from wgscovplot.features import EChartsFeature
from wgscovplot.tools.mosdepth import MosdepthDepthInfo


class TemplateDB(BaseModel):
    """Singleton wgscovplot template DB class containing all data to populate HTML plots

    The single instance of this class will be serialized into a JS object for populating plots in the final HTML output.
    """
    samples: List[str]
    ref_seq: str
    depths: Dict[str, str]
    amplicon_depths: Dict[str, List[Dict]]
    mosdepth_info: Dict[str, MosdepthDepthInfo]
    variants: Dict[str, List[Dict[str, str]]]
    show_amplicons: bool = False
    show_genes: bool = False
    segment_virus: bool = False
    low_coverage_threshold: int = 10
    max_depth: int = 10000
    echart_features: List[EChartsFeature]
