import logging
from pathlib import Path

from Bio import Entrez

import wgscovplot.flu
from wgscovplot import flu
from wgscovplot.db import DB, SegmentedGenomeDB
from wgscovplot.features import build_echarts_features_array
from wgscovplot.io import get_ref_seq_and_annotation, write_html_coverage_plot
from wgscovplot.rtpcr import flu_rtpcr_matches
from wgscovplot.tools import mosdepth, variants
from wgscovplot.tools.mosdepth.flu import get_flu_mosdepth_info

logger = logging.getLogger(__name__)
Entrez.email = "wgscovplot@github.com"


def run(
    input_dir: Path,
    low_coverage_threshold: int,
    primer_seq_path: Path,
    edit_distance: int,
    output_html: Path,
    compress_depths: bool = True,
) -> None:
    # Try to get sample-segment top references from "*.topsegments.csv" files
    # If no "*.topsegments.csv" files are found, then the virus is not segmented
    sample_top_references = flu.get_sample_top_references(input_dir)
    if sample_top_references:
        logger.info(f"Found segmented virus data in '{input_dir}'")
        db = build_db_segmented(
            input_dir=input_dir,
            sample_top_references=sample_top_references,
            primer_seq_path=primer_seq_path,
            low_coverage_threshold=low_coverage_threshold,
            edit_distance=edit_distance,
            compress_depths=compress_depths,
        )
    else:
        logger.info(f"Found unsegmented analysis data in '{input_dir}'")
        db = build_db(
            input_dir=input_dir, low_coverage_threshold=low_coverage_threshold, compress_depths=compress_depths
        )
    # Write coverage plot to HTML file
    write_html_coverage_plot(output_html=output_html, db=db)
    logger.info(f"Coverage plot written to '{output_html}'")


def build_db(
    input_dir: Path,
    low_coverage_threshold: int,
    compress_depths: bool = True,
) -> DB:
    ref_seq, gene_features = get_ref_seq_and_annotation(input_dir)
    # Get amplicon data
    amplicon_depths = mosdepth.get_amplicon_depths(input_dir)
    region_amplicon_data = mosdepth.get_region_amplicon(input_dir)
    if not (amplicon_depths and region_amplicon_data):
        logging.warning(f'No Mosdepth region BED file with amplicon data found in "{input_dir}"')
    echarts_features = build_echarts_features_array(gene_features, region_amplicon_data)
    # Get the list of samples name
    samples = mosdepth.get_samples_name(input_dir, is_genome_segmented=False)
    # Get coverage statistics information for all samples
    logger.info(f"Getting coverage depth info from {input_dir}")
    mosdepth_info, coverage_depths = mosdepth.get_info(input_dir, low_coverage_threshold=low_coverage_threshold)
    logger.info(f"{mosdepth_info=}")
    # Compress coverage depth arrays to base64 encoded strings if requested
    if compress_depths:
        coverage_depths_str = mosdepth.get_base64_encoded_depth_arrays(coverage_depths)
    else:
        coverage_depths_str = {sample: arr.tolist() for sample, arr in coverage_depths.items()}
    # Get Variant data
    samples_variants_info = variants.get_info(input_dir)
    variants_data = {
        sample: df_variants.to_dict(orient="records") for sample, df_variants in samples_variants_info.items()
    }
    return DB(
        samples=samples,
        ref_seq=ref_seq,
        depths=coverage_depths_str,
        amplicon_depths=amplicon_depths,
        mosdepth_info=mosdepth_info,
        variants=variants_data,
        low_coverage_threshold=low_coverage_threshold,
        echart_features=echarts_features,
    )


def build_db_segmented(
    input_dir: Path,
    sample_top_references: list[wgscovplot.flu.SampleSegmentRef],
    primer_seq_path: Path,
    low_coverage_threshold: int,
    edit_distance: int,
    compress_depths: bool = True,
) -> SegmentedGenomeDB:
    # Get list of samples name
    samples = mosdepth.get_samples_name(input_dir, is_genome_segmented=True)
    logger.info(f"Found {len(samples)} samples in '{input_dir}'")
    # Get reference id for each segment of each sample
    sample_segment_ref_id = flu.get_segments_ref_id(sample_top_references)
    logger.info(f"Found {len(sample_segment_ref_id)} segments in '{input_dir}'")
    # Get reference seq for each segment of each sample
    sample_segment_seq = flu.get_sample_segment_seqs(input_dir, sample_top_references)
    logger.info(f"Found {len(sample_segment_seq)} segment sequences in '{input_dir}'")
    # Get coverage depth for each segment of each sample
    mosdepth_info, coverage_depths = get_flu_mosdepth_info(input_dir, sample_top_references, low_coverage_threshold)
    # Compress coverage depth arrays to base64 encoded strings if requested
    if compress_depths:
        coverage_depths_str = mosdepth.get_base64_encoded_depth_arrays(coverage_depths)
    else:
        coverage_depths_str = {
            sample: {segment: arr.tolist() for segment, arr in segment_depths.items()}
            for sample, segment_depths in coverage_depths.items()
        }
    # Get variant info for each segment of each sample
    variants_data = variants.get_nf_flu_variant_info(
        basedir=input_dir, req_seq=sample_segment_seq, sample_top_references=sample_top_references
    )
    primer_matches: dict = {}
    if primer_seq_path is not None:
        # Get consensus sequence
        consensus_seq = flu.get_sample_segment_seqs(input_dir, sample_top_references)
        primer_matches = flu_rtpcr_matches(primer_seq_path, consensus_seq, edit_distance)
    return SegmentedGenomeDB(
        samples=samples,
        segments=flu.SEGMENTS,
        segments_ref_id=sample_segment_ref_id,
        segments_ref_seq=sample_segment_seq,
        depths=coverage_depths_str,
        variants=variants_data,
        mosdepth_info=mosdepth_info,
        primer_matches=primer_matches,
        low_coverage_threshold=low_coverage_threshold,
    )
