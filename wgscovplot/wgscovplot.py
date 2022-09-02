import logging
import markdown
from pathlib import Path

from Bio import SeqIO, Entrez
from wgscovplot.tools import variants, mosdepth
from wgscovplot.prepare_data import get_gene_amplicon_feature, write_html_coverage_plot, \
    write_html_coverage_plot_segment_virus, stat_info, get_primer_data

logger = logging.getLogger(__name__)
Entrez.email = "wgscovplot@github.com"


def run(input_dir: Path, ref_seq_path: Path, gene_feature_path: Path, ncbi_accession_id: str,
        low_coverage_threshold: int, amplicon: bool, gene_feature: bool, segment_virus: bool,
        gene_misc_feature: bool, primer_seq_path: Path, edit_distance: int, dev: bool, output_html: Path) -> None:
    # Read README.md
    dirpath = Path(__file__).parent.parent
    readme = dirpath / 'README.md'
    with open(readme, "r", encoding="utf-8") as input_file:
        text = input_file.read()
    about_html = markdown.markdown(text, extensions=['tables', 'nl2br', 'extra', 'md_in_html'])

    if segment_virus:
        # Get list of samples name
        samples_name = mosdepth.get_samples_name(input_dir, segment_virus)
        # Find files sample.topsegments.csv for each sample
        sample_top_references = mosdepth.get_sample_top_references(input_dir)
        # Get list of segments
        segments_name = mosdepth.get_segments_name(sample_top_references)
        # Get reference id for each segment of each sample
        ref_id = mosdepth.get_segments_ref_id(segments_name, sample_top_references)
        # Get reference seq for each segment of each sample
        ref_seq = mosdepth.get_segments_ref_seq(input_dir, segments_name, sample_top_references)
        # Get coverage depth for each segment of each sample
        depths_data = mosdepth.get_segments_depth(input_dir, segments_name, sample_top_references)
        # Get variant info for each segment of each sample
        variants_data = variants.get_segments_variants(input_dir, segments_name, sample_top_references)
        # Get summary of coverage statistics
        summary_info = mosdepth.get_flu_info(input_dir, samples_name, segments_name,
                                             sample_top_references, low_coverage_threshold)
        # Get low coverage regions
        low_coverage_regions = mosdepth.get_low_coverage_regions(input_dir, segments_name, sample_top_references,
                                                                 low_coverage_threshold)
        primer_data = {}
        if primer_seq_path is not None:
            primer_data = get_primer_data(primer_seq_path, edit_distance, ref_seq)
        write_html_coverage_plot_segment_virus(samples_name=samples_name,
                                               segments_name=segments_name,
                                               depths_data=depths_data,
                                               variants_data=variants_data,
                                               ref_seq=ref_seq,
                                               ref_id=ref_id,
                                               summary_info=summary_info,
                                               low_coverage_regions=low_coverage_regions,
                                               low_coverage_threshold=low_coverage_threshold,
                                               primer_data=primer_data,
                                               about_html=about_html,
                                               output_html=output_html)
    else:
        ref_name = mosdepth.get_refseq_name(input_dir)
        if ref_name != '' and ncbi_accession_id == "" and ref_seq_path is None:
            ncbi_accession_id = ref_name
        if ref_seq_path is None and ncbi_accession_id == "":
            logger.error('Please provide reference sequence --ref-seq /path/to/reference_sequence.fasta '
                         'OR provide correct NCBI Accession ID with option --ncbi-accession-id')
            exit(1)
        # Parse reference sequence
        elif ref_seq_path is not None:
            with open(ref_seq_path) as fh:
                for name, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
                    ref_seq = seq
        else:
            try:
                logger.info(f'Fetching reference sequence with accession_id "{ncbi_accession_id}" from NCBI database')
                with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text",
                                   id=ncbi_accession_id) as fasta_handle:
                    for name, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                        ref_seq = seq
            except Exception as e:
                logger.error(e, exc_info=True)
                logger.error(
                    f'Error! can not fetch "{ncbi_accession_id}" please correct accession id by providing option '
                    f'--ncbi-accession-id OR '
                    f'provide option --ref-seq /path/to/reference_sequence.fasta ')
                exit(1)
        # Get the list of samples name
        samples_name = mosdepth.get_samples_name(input_dir, segment_virus)

        # Get amplicon data
        if amplicon:
            region_amplicon_depth_data = mosdepth.get_depth_amplicon(input_dir)
            region_amplicon_data = mosdepth.get_region_amplicon(input_dir)
            if not (region_amplicon_depth_data and region_amplicon_data):
                logging.warning('No amplicon data found')
                amplicon = False
                region_amplicon_depth_data = {}
                region_amplicon_data = {}
        else:
            region_amplicon_depth_data = {}
            region_amplicon_data = {}

        # Get gene/amplicon feature
        if gene_feature and gene_feature_path is None and ncbi_accession_id == "":
            logger.error('If you want to plot gene features, '
                         'please provide genbank file for gene features, '
                         'option --genbank /path/to/genbank.gb '
                         'OR provide NCBI Accession ID with option --ncbi-accession-id')
            exit(1)
        gene_amplicon_feature_data = get_gene_amplicon_feature(gene_feature, gene_misc_feature,
                                                               gene_feature_path, ncbi_accession_id,
                                                               region_amplicon_data)

        # Get gene feature name, this is used for larger viral genome, use select2 to allow quick nanigation to feature
        # of interest
        gene_feature_name = []
        if len(gene_amplicon_feature_data):
            for feature in gene_amplicon_feature_data:
                gene_feature_name.append(feature["name"])

        # Get coverage statistics information for all samples
        mosdepth_info = mosdepth.get_info(input_dir, low_coverage_threshold=low_coverage_threshold)
        summary_info = stat_info(mosdepth_info, low_coverage_threshold=low_coverage_threshold)

        # Get Variant data
        samples_variants_info = variants.get_info(input_dir)
        variants_data = {}
        for sample, df_variants in samples_variants_info.items():
            variants_data[sample] = df_variants.to_dict(orient='records')

        # Get Depths data
        depths_data = mosdepth.get_depth(input_dir)

        # Get Coverage stat for summary information
        coverage_stat = {}
        for sample, coverage_info in mosdepth_info.items():
            coverage_stat[sample] = coverage_info.dict()

        # Write coverage plot to HTML file
        write_html_coverage_plot(samples_name=samples_name,
                                 depths_data=depths_data,
                                 variants_data=variants_data,
                                 ref_seq=ref_seq,
                                 summary_info=summary_info,
                                 gene_amplicon_feature_data=gene_amplicon_feature_data,
                                 gene_feature_name=gene_feature_name,
                                 about_html=about_html,
                                 output_html=output_html,
                                 region_amplicon_depth_data=region_amplicon_depth_data,
                                 amplicon=amplicon,
                                 gene_feature=gene_feature,
                                 low_coverage_threshold=low_coverage_threshold,
                                 dev=dev)
