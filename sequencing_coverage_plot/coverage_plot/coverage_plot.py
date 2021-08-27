import logging
import pandas as pd
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from enum import Enum
from Bio import SeqIO


class Resources(Enum):
    ECHARTS: str = "https://cdn.jsdelivr.net/npm/echarts@5.1.2/dist/echarts.min.js"
    GENE_FEATURE_COLOR =['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#FF33D3', '#b15928', '#0006fc']


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


def get_gene_feature(annotation: Path) -> list:
    
    gene_feature =[]
    for seq_record in SeqIO.parse(annotation, "genbank"):
        index = 0
        for seq_feature in seq_record.features:
            if (seq_feature.type != "CDS" and seq_feature.type != "source"):
                feature_dict = {
                    "name": "",
                    "value": [],
                    "itemStyle":{
                        "color":''
                    }
                }
                start_pos =  int(seq_feature.location.start)
                end_pos = int(seq_feature.location.end)
                strand = seq_feature.strand
                if seq_feature.type == "5'UTR" or seq_feature.type == "3'UTR":
                    feature_dict.update({"name": seq_feature.type})
                    feature_dict.update({"value":[index, start_pos, end_pos, 0, strand]})
                    feature_dict["itemStyle"].update({"color": Resources.GENE_FEATURE_COLOR.value[index]})
                else:

                    feature_name = seq_feature.qualifiers['gene'][0]
                    feature_dict.update({"name": feature_name})
                    feature_dict.update({"value":[index, start_pos, end_pos, 0, strand]})
                    feature_dict["itemStyle"].update({"color": Resources.GENE_FEATURE_COLOR.value[index]})
                index = index + 1
                gene_feature.append(feature_dict)
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
    template_file = render_env.get_template("stack_coverage_template.html")

    #coverage_stat = [['Sample 01', '4743.9X', '4743.9X', '99.6%', 114, 93, '1-10; 29800-29903', '1-2; 29813-29903'], ['Sample 02', '14743.9X', '5343.9X', '99.6%', 2514, 86, '1-10; 29800-29903', '1-2; 29813-29903; 1-2; 29813-29903; 1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903;1-2; 29813-29903']]
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(template_file.render(samples_name=samples_name,
                                        depth_data=depth_data,
                                        variant_data=variant_data,
                                        ref_seq=ref_seq,
                                        coverage_stat= coverage_stat,
                                        gene_feature=gene_feature,
                                        echarts_js=Resources.ECHARTS.value))


def prepare_data(samples_name: Path):
    df_samples = pd.read_table(samples_name, names=['coverage_depth_file', 'vcf_file'], index_col=0)
    depth_data = []
    variant_data = []
    coverage_stat =[]
    low = 10
    for sample in df_samples.index:
        logging.info(f'Preparing data for "{sample}"')
        df_coverage_depth = read_depths(df_samples.loc[sample, 'coverage_depth_file'])
        variant_info = {}
        df_vcf = parse_vcf(df_samples.loc[sample, 'vcf_file'])
        depth_data.append(df_coverage_depth.loc[:, 'depth'].to_list())
        for idx in df_vcf.index:
            variant_info[df_vcf.loc[idx, 'POS']] = [df_vcf.loc[idx, 'REF'], df_vcf.loc[idx, 'ALT']]
        variant_data.append(variant_info)

        ## Get Coverage Statistic for each samples ##
        low_depth = (df_coverage_depth.depth < 10)
        zero_depth = (df_coverage_depth.depth == 0)

        mean_cov = f'{df_coverage_depth.depth.mean():.1f}X'
        median_cov = f'{df_coverage_depth.depth.median():.1f}X'
        genome_cov = "{:.2%}".format((df_coverage_depth.depth >= low).sum() / df_coverage_depth.shape[0])
        pos_low_cov = low_depth.sum()
        pos_no_cov = zero_depth.sum()
        region_low_cov = get_interval_coords(df_coverage_depth, low-1)
        region_no_cov = get_interval_coords(df_coverage_depth, 0)

        coverage_stat.append([sample, mean_cov, median_cov, genome_cov, pos_low_cov, pos_no_cov, region_low_cov, region_no_cov])
    return df_samples.index.to_list(), depth_data, variant_data, coverage_stat
