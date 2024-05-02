from pydantic import BaseModel

from wgscovplot.features import EChartsFeature
from wgscovplot.tools.mosdepth import MosdepthDepthInfo
from wgscovplot.tools.mosdepth.flu import FluMosdepthDepthInfo


class DBMixin(BaseModel):
    samples: list[str]
    low_coverage_threshold: int = 10


class SegmentedGenomeDB(DBMixin):
    """Singleton class containing all data to populate HTML plots for segmented viruses.

    The single instance of this class will be serialized into a JS object for populating plots in the final HTML output.
    """

    segments: list[str]
    segments_ref_id: dict[str, dict[str, str]]
    segments_ref_seq: dict[str, dict[str, str]]
    depths: dict[str, dict[str, str]]
    mosdepth_info: dict[str, dict[str, FluMosdepthDepthInfo]]
    primer_matches: dict[str, dict[str, list[dict]]]
    variants: dict[str, dict[str, list[dict]]]


class DB(DBMixin):
    """Singleton class containing all data to populate HTML plots for viruses with an unsegmented genome.

    The single instance of this class will be serialized into a JS object for populating plots in the final HTML output.
    """

    ref_seq: str
    depths: dict[str, str]
    amplicon_depths: dict[str, list[dict]]
    mosdepth_info: dict[str, MosdepthDepthInfo]
    variants: dict[str, list[dict]]
    echart_features: list[EChartsFeature]
