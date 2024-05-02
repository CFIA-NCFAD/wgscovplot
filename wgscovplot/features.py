from itertools import cycle
from typing import Optional

from pydantic import BaseModel

from wgscovplot.colors import AmpliconColour, palette
from wgscovplot.config import FeaturesProps
from wgscovplot.util import overlap


class FeatureCoords(BaseModel):
    start: int = 0
    end: int = 0
    level: int = 0


class Feature(BaseModel):
    start: int
    end: int
    name: str
    strand: int = 1


class EChartsFeatureValue(BaseModel):
    idx: int = 0
    start: int = 1
    end: int = 1
    level: int = 0
    strand: int = 1
    rotate: float = 0.0
    type: str = "gene"


class EChartsItemStyle(BaseModel):
    color: str


class EChartsFeature(BaseModel):
    name: str
    value: EChartsFeatureValue
    itemStyle: EChartsItemStyle  # noqa: N815


def build_echarts_features_array(
    gene_features: Optional[list[Feature]],
    amplicon_features: Optional[list[Feature]],
    fp: Optional[FeaturesProps] = None,
) -> list[EChartsFeature]:
    """Build a list of gene/amplicon feature properties for ECharts

    Args:
        gene_features: List of gene Feature objects parsed from GFF or Genbank
        amplicon_features: List of amplicon Feature objects parsed from amplicon BED file
        fp: Feature properties for plotting with ECharts

    """
    if fp is None:
        fp = FeaturesProps()
    out: list[EChartsFeature] = []
    colour_cycle = cycle(palette)
    fcminus = FeatureCoords()
    fcplus = FeatureCoords()
    next_plus_strand_level = fp.plus_strand_level + fp.rec_items_height + fp.gene_feature_padding
    next_minus_strand_level = fp.minus_strand_level + fp.rec_items_height + fp.gene_feature_padding
    if gene_features:
        colour_cycle = cycle(colour_cycle)
        for feature in gene_features:
            start_pos: int = feature.start
            end_pos: int = feature.end
            strand = feature.strand
            if strand == 1:
                if overlap(fcplus.start, fcplus.end, start_pos, end_pos):
                    level = next_plus_strand_level
                    if fcplus.level == level:
                        level = fp.plus_strand_level
                else:
                    level = fp.plus_strand_level
                fcplus = FeatureCoords(start=start_pos, end=end_pos, level=level)
            else:
                if overlap(fcminus.start, fcminus.end, start_pos, end_pos):
                    level = next_minus_strand_level
                    if fcminus.level == level:
                        level = fp.minus_strand_level
                else:
                    level = fp.minus_strand_level
                fcminus = FeatureCoords(start=start_pos, end=end_pos, level=level)
            out.append(
                EChartsFeature(
                    name=feature.name,
                    value=EChartsFeatureValue(
                        idx=len(out),
                        start=start_pos,
                        end=end_pos,
                        level=level,
                        strand=strand,
                        rotate=0.5 if strand == 1 else -0.5,
                        type="gene",
                    ),
                    itemStyle=EChartsItemStyle(color=next(colour_cycle)),
                )
            )
    if amplicon_features:
        for feature in amplicon_features:
            primer_pool = int(feature.name.split("_")[-1])
            offset = fp.amplicon_offset if gene_features is not None else 0
            if primer_pool % 2 == 0:
                level = fp.amplicon_pool2_level + offset
                amplicon_color = AmpliconColour.pool2
            else:
                level = fp.amplicon_pool1_level + offset
                amplicon_color = AmpliconColour.pool1
            out.append(
                EChartsFeature(
                    name=feature.name,
                    value=EChartsFeatureValue(
                        idx=len(out),
                        start=feature.start,
                        end=feature.end,
                        level=level,
                        strand=feature.strand,
                        rotate=0,
                        type="amplicon",
                    ),
                    itemStyle=EChartsItemStyle(color=amplicon_color),
                )
            )
    return out
