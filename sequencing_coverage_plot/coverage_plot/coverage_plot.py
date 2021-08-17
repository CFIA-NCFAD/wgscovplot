import logging

import pandas as pd
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from enum import Enum


class Resources(Enum):
    ECHARTS: str = "https://cdn.jsdelivr.net/npm/echarts@5.1.2/dist/echarts.min.js"


def read_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['sample_name', 'reference', 'pos', 'depth'],
                       header=None)
    return df


def parse_vcf(vcf_path) -> pd.DataFrame:
    df = pd.read_table(vcf_path,
                       comment='#',
                       header=None,
                       names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample'])
    return df


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


def write_html_coverage_plot(samples_name: list,
                             depth_data: list,
                             variant_data: list,
                             ref_seq: str,
                             output_html: Path) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    template_file = render_env.get_template("stack_coverage_template.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(template_file.render(samples_name=samples_name,
                                        depth_data=depth_data,
                                        variant_data=variant_data,
                                        ref_seq=ref_seq,
                                        echarts_js=Resources.ECHARTS.value))


def prepare_data(samples_name: Path):
    df_samples = pd.read_table(samples_name, names=['coverage_depth_file', 'vcf_file'], index_col=0)
    depth_data = []
    variant_data = []
    for sample in df_samples.index:
        logging.info(f'Preparing data for "{sample}"')
        df_coverage_depth = read_depths(df_samples.loc[sample, 'coverage_depth_file'])
        variant_info = {}
        df_vcf = parse_vcf(df_samples.loc[sample, 'vcf_file'])
        depth_data.append(df_coverage_depth.loc[:, 'depth'].to_list())
        for idx in df_vcf.index:
            variant_info[df_vcf.loc[idx, 'POS']] = [df_vcf.loc[idx, 'REF'], df_vcf.loc[idx, 'ALT']]
        variant_data.append(variant_info)
    return df_samples.index.to_list(), depth_data, variant_data
