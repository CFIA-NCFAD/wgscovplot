import logging
import numpy as np
import sys
import requests
import pandas as pd
import os

from typing import Optional, Tuple
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO
from itertools import cycle

GENE_FEATURE_COLORS = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#fDBF6F', '#FF7F00',
                       '#CAB2D6',
                       '#6A3D9A', '#FF33D3', '#B15928', '#0006FC', '#2FB0EC', '#F3D742', '#2E9CE1', '#273D63',
                       '#980B92',
                       '#BBB873', '#EEC678', '#47E10B', '#E3139B', '#151179', '#293948', '#5F6005', '#FE24BE',
                       '#A7C36B',
                       '#D454DD', '#A68E2D', '#DB5AAC', '#405425', '#A608E4', '#533551', '#367521', '#64B875',
                       '#6DB011',
                       '#F5DD11', '#8A8517', '#F8E541', '#2D2A50', '#AAC3CC', '#C5D840', '#B79619', '#BBB2FE',
                       '#E37B03',
                       '#AFBB3E', '#74A110', '#E9877E', '#973F28', '#AFCA57', '#6E5EDE', '#B95FC3', '#C10AC8',
                       '#A59B67',
                       '#624F98', '#57A6AF', '#2650FB', '#94AAD1', '#5C1662', '#B8A1A1', '#104DB7', '#C6CBEE',
                       '#AA694D',
                       '#9B67DA', '#8DE7BC', '#866D49', '#72CEDC', '#574B7C', '#CD4474', '#593A60', '#2A6BB7',
                       '#286028',
                       '#6965EB', '#14CB29', '#956709', '#EB6D76', '#7A9A21', '#692C3C', '#AABBB5', '#1989AE',
                       '#D78DCC',
                       '#C43AAA', '#BBC863', '#E55F9D', '#741B13', '#6A7675', '#221A53', '#1804EC', '#D61D88',
                       '#1D50B3',
                       '#CF0E24', '#D791A9', '#0892FE', '#F5A865', '#91EBC2', '#9F650D', '#1B0A0F', '#1E9E88',
                       '#B42E38',
                       '#9710C9']

resources = {
    'echarts_js': 'https://cdn.jsdelivr.net/npm/echarts@5.1.2/dist/echarts.min.js',
    'jquery_js': 'https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js',
    'select2_css': 'https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css',
    'select2_js': 'https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/js/select2.min.js'
}


def read_regular_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['sample_name', 'reference', 'pos', 'depth'],
                       header=None)
    return df


def read_amplicon_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['reference', 'start', 'end', 'amplicon', 'depth'],
                       header=None)
    return df


def parse_vcf(vcf_file) -> Tuple[str, pd.DataFrame]:
    """Read VCF file into a DataFrame"""
    # https://github.com/peterk87/xlavir/blob/fbe11b4cef38bc291500c69d62b8599912c45887/xlavir/tools/variants.py#L211
    gzipped = vcf_file.endswith('.gz')
    with os.popen(f'zcat < {vcf_file.absolute()}') if gzipped else open(vcf_file) as fh:
        vcf_cols = []
        variant_caller = ''
        for line in fh:
            if line.startswith('##source='):
                variant_caller = line.strip().replace('##source=', '')
            if line.startswith('##nanopolish'):
                variant_caller = 'nanopolish'
            if line.startswith('#CHROM'):
                vcf_cols = line[1:].strip().split('\t')
                break
        df = pd.read_table(fh,
                           comment='#',
                           header=None,
                           names=vcf_cols)
        df = df[~df.duplicated(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER'], keep='first')]
    return variant_caller, df


def get_interval_coords(df, threshold=0):
    pos = df[df.depth <= threshold].pos
    coords = []
    for i, x in enumerate(pos):
        if coords:
            last = coords[-1][-1]
            if x == last + 1:
                coords[-1].append(x)
            else:
                coords.append([x])
        else:
            coords.append([x])
    return '; '.join([f'{xs[0]}-{xs[-1]}' for xs in coords])


def get_gene_feature(annotation: Path) -> list:
    gene_feature = []
    # number_of_colors = 100
    # color_pallet = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    # for i in range(number_of_colors)]
    # print (color_pallet)
    colour_cycle = cycle(GENE_FEATURE_COLORS)
    for seq_record in SeqIO.parse(annotation, "genbank"):
        index = 0
        for seq_feature in seq_record.features:
            if seq_feature.type in ["CDS", "source"]:
                continue
            if seq_feature.type in ["5'UTR", "3'UTR"]:
                feature_name = seq_feature.type
            else:
                feature_name = seq_feature.qualifiers['gene'][0]
            start_pos = int(seq_feature.location.start) + 1
            end_pos = int(seq_feature.location.end)
            strand = seq_feature.strand
            gene_feature.append(
                dict(name=feature_name,
                     value=[index, start_pos, end_pos, 0, strand],
                     itemStyle={"color": next(colour_cycle)})
            )
            index = index + 1
    return gene_feature


def write_html_coverage_plot(samples_name: list,
                             depth_data: list,
                             variant_data: list,
                             ref_seq: str,
                             coverage_stat: list,
                             gene_feature: list,
                             output_html: Path) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    template_file = render_env.get_template("covplot_template.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        logging.info('Retrieving JS and CSS for Coverage Plot')
        scripts_css = {}
        for k, v in resources.items():
            logging.info(f'Getting HTML resource "{k}" from "{v}"')
            scripts_css[k] = requests.get(v).text
        fout.write(template_file.render(samples_name=samples_name,
                                        depth_data=depth_data,
                                        variant_data=variant_data,
                                        coverage_stat=coverage_stat,
                                        gene_feature=gene_feature,
                                        ref_seq=ref_seq,
                                        ref_seq_length=len(ref_seq),
                                        **scripts_css))


def get_samples_name(df_samples) -> dict():
    samples_dict = {}
    for i, sample in enumerate(df_samples.index):
        logging.info(f'Preparing data for "{sample}"')
        samples_dict[i] = sample
    return samples_dict


def get_depth_data(df_samples, low=10):
    depth_data = {}
    coverage_stat = []
    for i, sample in enumerate(df_samples.index):
        df_coverage_depth = read_regular_depths(df_samples.loc[sample, 'coverage_depth_file'])
        depth_data[sample] = df_coverage_depth.loc[:, 'depth'].to_list()
        # Get Coverage Statistic for each samples
        low_depth = (df_coverage_depth.depth < 10)
        zero_depth = (df_coverage_depth.depth == 0)
        mean_cov = f'{df_coverage_depth.depth.mean():.1f}X'
        median_cov = f'{df_coverage_depth.depth.median():.1f}X'
        genome_cov = "{:.2%}".format((df_coverage_depth.depth >= low).sum() / df_coverage_depth.shape[0])
        pos_low_cov = low_depth.sum()
        pos_no_cov = zero_depth.sum()
        region_low_cov = get_interval_coords(df_coverage_depth, low - 1)
        region_no_cov = get_interval_coords(df_coverage_depth, 0)
        coverage_stat.append(
            [sample, mean_cov, median_cov, genome_cov, pos_low_cov, pos_no_cov, region_low_cov, region_no_cov])
    return depth_data, coverage_stat


def get_variant_data(df_samples) -> dict():
    variants_data = {}
    for i, sample in enumerate(df_samples.index):
        variant_info = {}
        if df_samples.loc[sample, 'vcf_file'] != 0:
            df_vcf = parse_vcf(df_samples.loc[sample, 'vcf_file'])[1]
            variant_info = {row.POS: (row.REF, row.ALT) for row in df_vcf.itertuples()}
        variants_data[sample] = variant_info
    return variants_data
