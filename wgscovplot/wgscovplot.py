import logging
from pathlib import Path

from Bio import Entrez

import wgscovplot.flu
import wgscovplot.tools.mosdepth.flu as flu
import wgscovplot.util as util
from wgscovplot.features import build_echarts_features_array
from wgscovplot.io import write_html_coverage_plot_segment_virus, write_html_coverage_plot, get_ref_seq_and_annotation, \
    TemplateHTML
from wgscovplot.db import TemplateDB
from wgscovplot.rtpcr import flu_rtpcr_matches
from wgscovplot.stats import cov_stats_to_html_table
from wgscovplot.tools import variants, mosdepth

logger = logging.getLogger(__name__)
Entrez.email = "wgscovplot@github.com"


def run(
        input_dir: Path,
        low_coverage_threshold: int,
        show_amplicon: bool,
        show_gene_features: bool,
        is_segmented: bool,
        primer_seq_path: Path,
        edit_distance: int,
        dev: bool,
        output_html: Path
) -> None:
    if is_segmented:
        # Get list of samples name
        samples = flu.get_samples_name(input_dir, is_segmented)
        # Find files sample.topsegments.csv for each sample
        sample_top_references = wgscovplot.flu.get_sample_top_references(input_dir)
        # Get list of segments
        segments = flu.get_segments(sample_top_references)
        # Get reference id for each segment of each sample
        ref_id = flu.get_segments_ref_id(sample_top_references)
        # Get reference seq for each segment of each sample
        ref_seq = flu.get_segments_ref_seq(input_dir, sample_top_references)
        # Get coverage depth for each segment of each sample
        depths = flu.get_segments_depth(input_dir, sample_top_references)
        # Get variant info for each segment of each sample
        variants_data = variants.get_segments_variants(
            input_dir,
            sample_top_references
        )
        # Get summary of coverage statistics
        cov_html_table = flu.get_flu_info(
            basedir=input_dir,
            samples=samples,
            segments=segments,
            sample_top_references=sample_top_references,
            low_coverage_threshold=low_coverage_threshold
        )
        # Get low coverage regions
        low_coverage_regions = flu.get_low_coverage_regions(input_dir, segments, sample_top_references,
                                                            low_coverage_threshold)
        primer_data = {}
        if primer_seq_path is not None:
            # TODO: sample consensus sequences need to be matched against NOT ref seqs!!!
            primer_data = flu_rtpcr_matches(primer_seq_path, ref_seq, edit_distance)
        write_html_coverage_plot_segment_virus(samples_name=samples,
                                               segments_name=segments,
                                               depths_data=depths,
                                               variants_data=variants_data,
                                               ref_seq=ref_seq,
                                               ref_id=ref_id,
                                               summary_info=cov_html_table,
                                               low_coverage_regions=low_coverage_regions,
                                               low_coverage_threshold=low_coverage_threshold,
                                               primer_data=primer_data,
                                               about_html=util.readme_to_html(),
                                               output_html=output_html)
    else:
        ref_seq, gene_features = get_ref_seq_and_annotation(input_dir, get_gene_features=show_gene_features)
        # Get amplicon data
        if show_amplicon:
            amplicon_depths = mosdepth.get_amplicon_depths(input_dir)
            region_amplicon_data = mosdepth.get_region_amplicon(input_dir)
            if not (amplicon_depths and region_amplicon_data):
                logging.warning(f'No Mosdepth region BED file with amplicon data found in "{input_dir}"')
                show_amplicon = False
                amplicon_depths = {}
                region_amplicon_data = {}
        else:
            amplicon_depths = {}
            region_amplicon_data = {}
        echarts_features = build_echarts_features_array(
            gene_features,
            region_amplicon_data
        )

        # Get the list of samples name
        samples = wgscovplot.tools.mosdepth.flu.get_samples_name(input_dir, is_segmented)

        # Get coverage statistics information for all samples
        mosdepth_info, sample_depths = mosdepth.get_info(input_dir, low_coverage_threshold=low_coverage_threshold)
        cov_html_table = cov_stats_to_html_table(mosdepth_info)

        # Get Variant data
        samples_variants_info = variants.get_info(input_dir)
        variants_data = {}
        for sample, df_variants in samples_variants_info.items():
            variants_data[sample] = df_variants.to_dict(orient='records')

        # encode sample coverage depth arrays in base64 for better compression vs dumping a regular list to JSON
        depths = mosdepth.get_base64_encoded_depth_arrays(sample_depths)
        # depths = {sample: ds.astype(int).tolist() for sample, ds in sample_depths.items()}

        db = TemplateDB(
            samples=samples,
            ref_seq=ref_seq,
            depths=depths,
            amplicon_depths=amplicon_depths,
            mosdepth_info=mosdepth_info,
            variants=variants_data,
            show_amplicon=show_amplicon,
            show_gene_features=show_gene_features,
            low_coverage_threshold=low_coverage_threshold,
            echart_features=echarts_features
        )
        html = TemplateHTML(
            about_html=util.readme_to_html(),
            cov_stats_html=cov_html_table,
        )
        # Write coverage plot to HTML file
        write_html_coverage_plot(
            output_html=output_html,
            db=db,
            html=html,
            dev=dev,
        )
