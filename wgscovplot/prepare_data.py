import logging
import math
import pandas as pd
from typing import Dict, List, Any
from pathlib import Path
from itertools import product
import edlib

from Bio.SeqFeature import SeqFeature
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO, Entrez
from itertools import cycle
from .resources import gene_feature_properties
from .colors import color_pallete, AmpliconColour
from wgscovplot.tools import mosdepth

logger = logging.getLogger(__name__)

NT_MAP = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}


def expand_degenerate_bases(seq: str) -> List[str]:
    return list(map("".join, product(*map(NT_MAP.get, seq))))


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
               'Mean Coverage Depth (X)',
               'Median Coverage Depth (X)', 'Ref Sequence Length (bp)']
    df.columns = headers
    df.drop(columns=['Low Coverage Threshold'], inplace=True)
    for index in df.index:
        for col in df.columns:
            df.loc[index, col] = "{:.2%}".format(df.loc[index, col][1]) if col == '% Genome Coverage >= ' + str(
                low_coverage_threshold) + 'X' else \
                df.loc[index, col][1]
    df.sort_values(by=['Sample'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.index = df.index + 1
    return df.to_html(classes="table table-striped table-hover table-bordered table-responsive-md",
                      float_format=lambda x: f'{x:0.2f}', justify="left", table_id="summary-coverage-stat")


def get_gene_amplicon_feature(gene_feature: bool, gene_misc_feature: bool,
                              annotation: Path, ncbi_accession_id: str,
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
            except Exception as e:
                logger.error(e, exc_info=True)
                logger.error(
                    f'Error! can not fetch "{ncbi_accession_id}" please correct accession id '
                    f'by providing option --ncbi-accession-id OR '
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
                             "rotate": 0.5 if strand == 1 else -0.5,
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
                         "rotate": 0,
                         "type": "amplicon_feature"
                     },
                     itemStyle={"color": amplicon_color})
            )
    return feature_data


def get_primer_data(primer_seq_path: Path, edit_distance: int, ref_seq: Dict[str, Dict[str, str]]) -> \
        Dict[str, Dict[str, Dict]]:
    primer_data = {}
    primer_seq_record = []
    for record in SeqIO.parse(open(primer_seq_path), 'fasta'):
        primer_seq_record.append([record.id, record.seq])
    for sample_key in ref_seq.keys():
        primer_data[sample_key] = {}
        for segment_key in ref_seq[sample_key].keys():
            primer_data[sample_key][segment_key] = []
            if ref_seq[sample_key][segment_key] != '':
                for idx, record in enumerate(primer_seq_record):
                    seqs = expand_degenerate_bases(record[1])
                    align_result = []
                    for seq in seqs:
                        align = edlib.align(seq, ref_seq[sample_key][segment_key],
                                            mode="HW", task="path", k=edit_distance)
                        align['query_seq'] = seq
                        if align['editDistance'] != -1:  # keep sequence has alignment result
                            align_result.append(align)
                    if len(align_result):
                        align_min_ed = min(align_result, key=lambda x: x['editDistance'])
                        other_locations = []
                        for location in align_min_ed['locations'][1:]:
                            loc_start = location[0]
                            loc_end = location[1]
                            # covert to 1-based index
                            if loc_start is not None:
                                loc_start = loc_start + 1
                            if loc_end is not None:
                                loc_end = loc_end + 1
                            other_locations.append((loc_start, loc_end))
                        viz_result = edlib.getNiceAlignment(align_min_ed, align_min_ed['query_seq'],
                                                            ref_seq[sample_key][segment_key])
                        primer_data[sample_key][segment_key].append(dict(name=record[0],
                                                                         query_aligned=str(viz_result['query_aligned']),
                                                                         target_aligned=viz_result['target_aligned'],
                                                                         matched_aligned=viz_result['matched_aligned'],
                                                                         # nice_viz="\n".join(viz_result.values()),
                                                                         cigar=align_min_ed['cigar'],
                                                                         edit_distance=align_min_ed['editDistance'],
                                                                         start=align_min_ed['locations'][0][0],
                                                                         end=align_min_ed['locations'][0][1],
                                                                         other_locations=", ".join(
                                                                             map(str, other_locations))))
    return primer_data


def write_html_coverage_plot(samples_name: List[str],
                             depths_data: Dict[str, List],
                             variants_data: Dict[str, List],
                             ref_seq: str,
                             summary_info: str,
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
                                        summary_info=summary_info,
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
                                           variants_data: Dict[str, Dict[str, Dict]],
                                           ref_seq: Dict[str, Dict[str, str]],
                                           ref_id: Dict[str, Dict[str, str]],
                                           summary_info: str,
                                           low_coverage_regions: Dict[str, Dict[str, str]],
                                           low_coverage_threshold: int,
                                           primer_data: Dict[str, Dict[str, Dict]],
                                           about_html: str,
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
                                        variants_data=variants_data,
                                        ref_seq=ref_seq,
                                        ref_id=ref_id,
                                        summary_info=summary_info,
                                        low_coverage_regions=low_coverage_regions,
                                        low_coverage_threshold=low_coverage_threshold,
                                        primer_data=primer_data,
                                        about_html=about_html))
