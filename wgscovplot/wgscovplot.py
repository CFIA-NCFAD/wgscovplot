import logging
import requests
import math
import markdown
import pandas as pd
from typing import Dict, List, Any
from pathlib import Path

from Bio.SeqFeature import SeqFeature
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO, Entrez
from itertools import cycle
from .resources import cdn_resources, gene_feature_properties
from .colors import color_pallete, AmpliconColour
from wgscovplot.tools import variants, mosdepth

logger = logging.getLogger(__name__)
Entrez.email = "nhhaidee@gmail.com"


def run(input_dir: Path, ref_seq: Path, genbank: Path, ncbi_accession_id: str, amplicon: bool, gene_feature: bool,
        gene_misc_feature: bool, output_html: Path) -> None:
    if ref_seq is None and ncbi_accession_id == "":
        logger.error('Please provide reference sequence --ref-seq /path/to/reference_sequence.fasta '
                     'OR provide NCBI Accession ID with option --ncbi-accession-id')
        exit(1)
    # Parse reference sequence
    elif ref_seq is not None:
        with open(ref_seq) as fh:
            for name, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
                ref_seq = seq
    else:
        try:
            with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=ncbi_accession_id) as fasta_handle:
                for name, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                    ref_seq = seq
        except:
            logger.error(f'Error! can not fetch "{ncbi_accession_id}" please correct accession id OR '
                         f'provide option --ref-seq /path/to/reference_sequence.fasta ')
            exit(1)
    # Get the list of samples name
    samples_name = mosdepth.get_samples_name(input_dir)

    # Get amplicon data
    if amplicon:
        amplicon_depths_data = mosdepth.get_depth_amplicon(input_dir)
        amplicon_regions_data = mosdepth.get_region_amplicon(input_dir)
        if not (amplicon_depths_data and amplicon_regions_data):
            logging.warning('No amplicon data found')
            amplicon = False
    else:
        amplicon_depths_data = {}
        amplicon_regions_data = {}

    # Get gene/amplicon feature
    if gene_feature and genbank is None and ncbi_accession_id == "":
        logger.error('If you want to plot gene features, please provide genbank file for gene features, '
                     'option --genbank /path/to/genbank.gb OR provide NCBI Accession ID with option --ncbi-accession-id')
        exit(1)
    gene_feature_data = get_gene_feature(gene_feature, gene_misc_feature, genbank, ncbi_accession_id, amplicon_regions_data)

    # Get coverage statistics information for all samples
    mosdepth_info = mosdepth.get_info(input_dir, low_coverage_threshold=10)
    sample_stat_info = stat_info(mosdepth_info)

    samples_variants_info = variants.get_info(input_dir)

    # Get Variant matrix using for Variant Heatmap
    df_variants = variants.to_dataframe(samples_variants_info.values())
    mutation = []
    variant_matrix_data = []
    if 'Mutation' in df_variants.columns:
        df_varmap = variants.to_variant_pivot_table(df_variants)
        for i, sample in enumerate(samples_name):
            for j, mutation_name in enumerate(df_varmap.columns):
                if sample in df_varmap.index:
                    variant_matrix_data.append([j, i, df_varmap.loc[sample, mutation_name]])
                else:
                    variant_matrix_data.append([j, i, 0.0])
        mutation = df_varmap.columns.tolist()
    # Get Variant data
    variants_data = {}
    for sample, df_variants in samples_variants_info.items():
        variants_data[sample] = df_variants.to_dict(orient='records')

    # Get Depths data
    depths_data = mosdepth.get_depth(input_dir)

    # Get Coverage stat for summary inforamtion
    coverage_stat = {}
    for sample, coverage_info in mosdepth_info.items():
        coverage_stat[sample] = coverage_info.dict()

    # Read README.md
    dirpath = Path(__file__).parent
    readme = dirpath / 'readme/README.md'
    with open(readme, "r", encoding="utf-8") as input_file:
        text = input_file.read()
    about_html = markdown.markdown(text, extensions=['tables', 'nl2br', 'extra', 'md_in_html'])

    # Write coverage plot to HTML file
    write_html_coverage_plot(samples_name=samples_name,
                             depth_data=depths_data,
                             variant_data=variants_data,
                             ref_seq=ref_seq,
                             coverage_stat=sample_stat_info,
                             gene_feature_data=gene_feature_data,
                             gene_feature=gene_feature,
                             amplicon_data=amplicon_depths_data,
                             variant_matrix=variant_matrix_data,
                             mutation=mutation,
                             amplicon=amplicon,
                             about_html=about_html,
                             output_html=output_html)


def stat_info(sample_depth_info: Dict[str, mosdepth.MosdepthDepthInfo]) -> str:
    df = pd.DataFrame(sample_depth_info.values())
    headers = ['Sample', '# 0 Coverage Positions', '0 Coverage Regions', 'Low Coverage Threshold',
               '# Low Coverage Positions (< 10X)', 'Low Coverage Regions (< 10X)', '% Genome Coverage >= 10X',
               'Mean Coverage Depth',
               'Median Coverage Depth', 'Ref Sequence Length']
    df.columns = headers
    df.drop(columns=['Low Coverage Threshold'], inplace=True)
    for index in df.index:
        for col in df.columns:
            df.loc[index, col] = "{:.2%}".format(df.loc[index, col][1]) if col == '% Genome Coverage >= 10X' else \
                df.loc[index, col][1]
    df.sort_values(by=['Sample'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df.to_html(classes="table table-striped table-hover table-bordered table-responsive-md",
                      float_format=lambda x: f'{x:0.2f}')


def overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 <= start2 <= end1 or start1 <= end2 <= end1


def max_depth(depth_data: Dict[str, List]) -> int:
    max_value = 0
    for key, values in depth_data.items():
        if max_value <= max(values):
            max_value = max(values)
    return math.ceil(max_value * 1.5)


def get_gene_feature(gene_feature: bool, gene_misc_feature: bool, annotation: Path, ncbi_accession_id: str,
                     amplicon_regions: Dict[str, List]) -> List[Dict[str, Any]]:
    gene_feature_data = []
    if gene_feature:
        colour_cycle = cycle(color_pallete)
        minus_strains_list = [0, 0, 0]
        plus_strains_list = [0, 0, 0]
        if annotation is not None:
            genbank_handle = annotation
        else:
            try:
                genbank_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_accession_id)
            except:
                logger.error(f'Error! can not fetch "{ncbi_accession_id}" please correct accession id OR '
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
                gene_feature_data.append(
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
    if amplicon_regions:
        for amplicon_name, region in amplicon_regions.items():
            index = int(amplicon_name.split('_')[-1])
            gene_feature_len = len(gene_feature_data)
            if index % 2:
                level = 95 if gene_feature else 15
                amplicon_color = AmpliconColour.pool2.value
            else:
                level = 80 if gene_feature else 0
                amplicon_color = AmpliconColour.pool1.value
            gene_feature_data.append(
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
    return gene_feature_data


def write_html_coverage_plot(samples_name: List[str],
                             depth_data: Dict[str, List],
                             variant_data: Dict[str, List],
                             ref_seq: str,
                             coverage_stat: str,
                             gene_feature_data: List[Dict[str, Any]],
                             about_html: str,
                             output_html: Path,
                             amplicon_data: Dict[str, List],
                             variant_matrix: List[List],
                             mutation: List[str],
                             amplicon: bool = False,
                             gene_feature: bool = False
                             ) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    template_file = render_env.get_template("wgscovplot_template.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        logging.info('Retrieving JS and CSS resources for Coverage Plot')
        scripts_css = {}
        for k, v in cdn_resources.items():
            logging.info(f'Getting HTML resource "{k}" from "{v}"')
            scripts_css[k] = requests.get(v).text
        fout.write(template_file.render(gene_feature_data=gene_feature_data,
                                        gene_feature=gene_feature,
                                        amplicon_data=amplicon_data,
                                        amplicon=amplicon,
                                        samples_name=samples_name,
                                        depth_data=depth_data,
                                        variant_data=variant_data,
                                        coverage_stat=coverage_stat,
                                        ref_seq=ref_seq,
                                        ref_seq_length=len(ref_seq),
                                        variant_matrix=variant_matrix,
                                        mutation=mutation,
                                        about_html=about_html,
                                        max_depth=max_depth(depth_data),
                                        **scripts_css))
