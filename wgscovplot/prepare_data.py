import logging
import math
import pandas as pd
from typing import Dict, List, Any
from pathlib import Path

from Bio.SeqFeature import SeqFeature
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO, Entrez
from itertools import cycle
from .resources import gene_feature_properties
from .colors import color_pallete, AmpliconColour
from wgscovplot.tools import mosdepth

logger = logging.getLogger(__name__)


def overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 <= start2 <= end1 or start1 <= end2 <= end1


def max_depth(depth_data: Dict[str, List]) -> int:
    max_value = 0
    for key, values in depth_data.items():
        if max_value <= max(values):
            max_value = max(values)
    return math.ceil(max_value * 1.5)


def stat_info(sample_depth_info: Dict[str, mosdepth.MosdepthDepthInfo], low_coverage_threshold: int) -> str:
    df = pd.DataFrame(sample_depth_info.values())
    headers = ['Sample', '# 0 Coverage Positions', '0 Coverage Regions', 'Low Coverage Threshold',
               '# Low Coverage Positions (< ' + str(low_coverage_threshold) + 'X)',
               'Low Coverage Regions (< ' + str(low_coverage_threshold) + 'X)',
               '% Genome Coverage >= ' + str(low_coverage_threshold) + 'X',
               'Mean Coverage Depth',
               'Median Coverage Depth', 'Ref Sequence Length']
    df.columns = headers
    df.drop(columns=['Low Coverage Threshold'], inplace=True)
    for index in df.index:
        for col in df.columns:
            df.loc[index, col] = "{:.2%}".format(df.loc[index, col][1]) if col == '% Genome Coverage >= ' + str(
                low_coverage_threshold) + 'X' else \
                df.loc[index, col][1]
    df.sort_values(by=['Sample'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df.to_html(classes="table table-striped table-hover table-bordered table-responsive-md",
                      float_format=lambda x: f'{x:0.2f}')


def get_gene_amplicon_feature(gene_feature: bool, gene_misc_feature: bool, annotation: Path, ncbi_accession_id: str,
                              region_amplicon_data: Dict[str, List]) -> List[Dict[str, Any]]:
    feature_data = []
    if gene_feature:
        colour_cycle = cycle(color_pallete)
        minus_strains_list = [0, 0, 0]
        plus_strains_list = [0, 0, 0]
        if annotation is not None:
            genbank_handle = annotation
        else:
            try:
                logger.info(f'Fetching Genbank file with accession_id "{ncbi_accession_id}" from NCBI database')
                genbank_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_accession_id)
            except:
                logger.error(
                    f'Error! can not fetch "{ncbi_accession_id}" please correct accession id by providing option --ncbi-accession-id OR '
                    f'provide option --genbank /path/to/genbank.gb ')

                exit(1)
        if gene_misc_feature:
            skip_feature = ["CDS", "source", "repeat_region"]
        else:
            skip_feature = ["CDS", "source", "repeat_region", "misc_feature"]
        for seq_record in SeqIO.parse(genbank_handle, "gb"):
            index = 0  # the index must be continuous for data handling with Echarts
            seq_feature: SeqFeature
            for seq_feature in seq_record.features:
                if seq_feature.type in skip_feature:
                    continue
                if seq_feature.type in ["5'UTR", "3'UTR"]:
                    feature_name = seq_feature.type
                else:
                    if seq_feature.qualifiers.get('gene'):
                        feature_name = seq_feature.qualifiers['gene'][0]
                    elif seq_feature.qualifiers.get('locus_tag'):
                        feature_name = seq_feature.qualifiers['locus_tag'][0]
                    else:
                        feature_name = seq_feature.id
                start_pos = int(seq_feature.location.start) + 1
                end_pos = int(seq_feature.location.end)
                strand = seq_feature.strand
                if strand == 1:
                    if overlap(plus_strains_list[0], plus_strains_list[1], start_pos, end_pos):
                        level = gene_feature_properties['plus_strand_level'] + gene_feature_properties[
                            'rec_items_height'] + 3
                        if plus_strains_list[2] == gene_feature_properties['plus_strand_level'] + \
                                gene_feature_properties['rec_items_height'] + 3:
                            level = gene_feature_properties['plus_strand_level']
                    else:
                        level = gene_feature_properties['plus_strand_level']
                    plus_strains_list = [start_pos, end_pos, level]
                else:
                    if overlap(minus_strains_list[0], minus_strains_list[1], start_pos, end_pos):
                        level = gene_feature_properties['minus_strand_level'] + gene_feature_properties[
                            'rec_items_height'] + 3
                        if minus_strains_list[2] == gene_feature_properties['minus_strand_level'] + \
                                gene_feature_properties['rec_items_height'] + 3:
                            level = gene_feature_properties['minus_strand_level']
                    else:
                        level = gene_feature_properties['minus_strand_level']
                    minus_strains_list = [start_pos, end_pos, level]
                feature_data.append(
                    dict(name=feature_name,
                         value={
                             "idx": index,
                             "start": start_pos,
                             "end": end_pos,
                             "level": level,
                             "strand": strand,
                             "type": "gene_feature"
                         },
                         itemStyle={"color": next(colour_cycle)})
                )
                index += 1
    if region_amplicon_data:
        for amplicon_name, region in region_amplicon_data.items():
            index = int(amplicon_name.split('_')[-1])
            gene_feature_len = len(feature_data)
            if index % 2:
                level = 95 if gene_feature else 15
                amplicon_color = AmpliconColour.pool2.value
            else:
                level = 80 if gene_feature else 0
                amplicon_color = AmpliconColour.pool1.value
            feature_data.append(
                dict(name=amplicon_name,
                     value={
                         "idx": gene_feature_len,
                         "start": region[0],
                         "end": region[1],
                         "level": level,
                         "strand": '',
                         "type": "amplicon_feature"
                     },
                     itemStyle={"color": amplicon_color})
            )
    return feature_data


def write_html_coverage_plot(samples_name: List[str],
                             depths_data: Dict[str, List],
                             variants_data: Dict[str, List],
                             ref_seq: str,
                             coverage_stat: str,
                             gene_amplicon_feature_data: List[Dict[str, Any]],
                             gene_feature_name: List[str],
                             about_html: str,
                             output_html: Path,
                             region_amplicon_depth_data: Dict[str, List],
                             amplicon: bool = False,
                             gene_feature: bool = False,
                             low_coverage_threshold: int = 10,
                             dev: bool = False
                             ) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    if dev:
        template_file = render_env.get_template("wgscovplot_dev_template.html")
    else:
        template_file = render_env.get_template("wgscovplot_prod_template.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(template_file.render(samples_name=samples_name,
                                        depths_data=depths_data,
                                        variants_data=variants_data,
                                        ref_seq=ref_seq,
                                        coverage_stat=coverage_stat,
                                        gene_amplicon_feature_data=gene_amplicon_feature_data,
                                        gene_feature_name=gene_feature_name,
                                        about_html=about_html,
                                        region_amplicon_depth_data=region_amplicon_depth_data,
                                        amplicon=amplicon,
                                        gene_feature=gene_feature,
                                        ref_seq_length=len(ref_seq),
                                        max_depth=max_depth(depths_data),
                                        low_coverage_threshold=low_coverage_threshold))


def write_html_coverage_plot_segment_virus(samples_name: List[str],
                                           segments_name: List[str],
                                           depths_data: Dict[str, Dict[str, List]],
                                           ref_seq: Dict[str, Dict[str, str]],
                                           output_html: Path,
                                           ) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    template_file = render_env.get_template("wgscovplot_flu_template.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(template_file.render(samples_name=samples_name,
                                        segments_name=segments_name,
                                        depths_data=depths_data,
                                        ref_seq=ref_seq))
