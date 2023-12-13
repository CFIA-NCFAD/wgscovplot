import logging
from pathlib import Path

from Bio import Entrez
import json

import wgscovplot.flu
import wgscovplot.tools.mosdepth.flu
from wgscovplot.features import build_echarts_features_array
from wgscovplot.io import write_html_coverage_plot, get_ref_seq_and_annotation
from wgscovplot.db import NonSegmentedGenomeDB, SegmentedGenomeDB
from wgscovplot.rtpcr import flu_rtpcr_matches
from wgscovplot.tools import variants, mosdepth

logger = logging.getLogger(__name__)
Entrez.email = "wgscovplot@github.com"


def run(
        input_dir: Path,
        low_coverage_threshold: int,
        primer_seq_path: Path,
        edit_distance: int,
        output_html: Path,
) -> None:
    sample_top_references = wgscovplot.flu.get_sample_top_references(input_dir)
    if sample_top_references:
        # Get list of samples name
        samples = wgscovplot.tools.mosdepth.get_samples_name(input_dir, is_genome_segmented=True)
        # Get list of segments
        segments = wgscovplot.flu.get_segments(sample_top_references)
        # Get reference id for each segment of each sample
        ref_id = wgscovplot.flu.get_segments_ref_id(segments, sample_top_references)
        # Get reference seq for each segment of each sample
        ref_seq = wgscovplot.flu.get_segments_ref_seq(input_dir, segments, sample_top_references)
        # Get coverage depth for each segment of each sample
        mosdepth_info, coverage_depths = wgscovplot.tools.mosdepth.flu.get_flu_mosdepth_info(input_dir, segments,
                                                                                             sample_top_references,
                                                                                             low_coverage_threshold)
        # Get variant info for each segment of each sample
        variants_data = variants.get_segments_variants(
            basedir=input_dir,
            segments=segments,
            req_seq=ref_seq,
            sample_top_references=sample_top_references
        )
        primer_matches = {}
        if primer_seq_path is not None:
            # Get consensus sequence
            consensus_seq = wgscovplot.flu.get_segments_consensus_seq(input_dir, segments, sample_top_references)
            primer_matches = flu_rtpcr_matches(primer_seq_path, consensus_seq, edit_distance)
        db = SegmentedGenomeDB(
            samples=samples,
            segments=segments,
            segments_ref_id=ref_id,
            segments_ref_seq=ref_seq,
            depths=coverage_depths,
            variants=variants_data,
            mosdepth_info=mosdepth_info,
            primer_matches=primer_matches,
            low_coverage_threshold=low_coverage_threshold,
        )

        write_html_coverage_plot(
            output_html=output_html,
            db=db
        )
    else:
        ref_seq, gene_features = get_ref_seq_and_annotation(input_dir)
        # Get amplicon data
        amplicon_depths = mosdepth.get_amplicon_depths(input_dir)
        region_amplicon_data = mosdepth.get_region_amplicon(input_dir)
        if not (amplicon_depths and region_amplicon_data):
            logging.warning(f'No Mosdepth region BED file with amplicon data found in "{input_dir}"')
        echarts_features = build_echarts_features_array(gene_features, region_amplicon_data)

        # Get the list of samples name
        samples = wgscovplot.tools.mosdepth.get_samples_name(input_dir, is_genome_segmented=False)

        # Get coverage statistics information for all samples
        mosdepth_info, coverage_depths = mosdepth.get_info(input_dir, low_coverage_threshold=low_coverage_threshold)

        # Get Variant data
        samples_variants_info = variants.get_info(input_dir)
        variants_data = {}
        for sample, df_variants in samples_variants_info.items():
            variants_data[sample] = df_variants.to_dict(orient='records')

        db = NonSegmentedGenomeDB(
            samples=samples,
            ref_seq=ref_seq,
            depths=coverage_depths,
            amplicon_depths=amplicon_depths,
            mosdepth_info=mosdepth_info,
            variants=variants_data,
            low_coverage_threshold=low_coverage_threshold,
            echart_features=echarts_features,
        )
        # Write coverage plot to HTML file
        write_html_coverage_plot(
            output_html=output_html,
            db=db
        )
