import logging
import markdown
from pathlib import Path

from Bio import SeqIO, Entrez
from wgscovplot.tools import variants, mosdepth
from wgscovplot.prepare_data import get_gene_amplicon_feature, write_html_coverage_plot, stat_info

logger = logging.getLogger(__name__)
Entrez.email = "wgscovplot@github.com"


def run(input_dir: Path, ref_seq: Path, genbank: Path, ncbi_accession_id: str, low_coverage_threshold: int,
        amplicon: bool, gene_feature: bool,
        gene_misc_feature: bool, output_html: Path) -> None:
    ref_name = mosdepth.get_refseq_name(input_dir)
    if ref_name != '' and ncbi_accession_id == "" and ref_seq is None:
        ncbi_accession_id = ref_name
    if ref_seq is None and ncbi_accession_id == "":
        logger.error('Please provide reference sequence --ref-seq /path/to/reference_sequence.fasta '
                     'OR provide correct NCBI Accession ID with option --ncbi-accession-id')
        exit(1)
    # Parse reference sequence
    elif ref_seq is not None:
        with open(ref_seq) as fh:
            for name, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
                ref_seq = seq
    else:
        try:
            logger.info(f'Fetching reference sequence with accession_id "{ncbi_accession_id}" from NCBI database')
            with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=ncbi_accession_id) as fasta_handle:
                for name, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                    ref_seq = seq
        except:
            logger.error(
                f'Error! can not fetch "{ncbi_accession_id}" please correct accession id by provding option --ncbi-accession-id OR '
                f'provide option --ref-seq /path/to/reference_sequence.fasta ')
            exit(1)
    # Get the list of samples name
    samples_name = mosdepth.get_samples_name(input_dir)

    # Get amplicon data
    if amplicon:
        region_amplicon_depth_data = mosdepth.get_depth_amplicon(input_dir)
        region_amplicon_data = mosdepth.get_region_amplicon(input_dir)
        if not (region_amplicon_depth_data and region_amplicon_data):
            logging.warning('No amplicon data found')
            amplicon = False
    else:
        region_amplicon_depth_data = {}
        region_amplicon_data = {}

    # Get gene/amplicon feature
    if gene_feature and genbank is None and ncbi_accession_id == "":
        logger.error('If you want to plot gene features, please provide genbank file for gene features, '
                     'option --genbank /path/to/genbank.gb OR provide NCBI Accession ID with option --ncbi-accession-id')
        exit(1)
    gene_amplicon_feature_data = get_gene_amplicon_feature(gene_feature, gene_misc_feature,
                                                           genbank, ncbi_accession_id, region_amplicon_data)

    # Get gene feature name, this is used for larger viral genome, use select2 to allow quick nanigation to feature
    # of interest
    gene_feature_name = []
    if len(gene_amplicon_feature_data):
        for feature in gene_amplicon_feature_data:
            gene_feature_name.append(feature["name"])

    # Get coverage statistics information for all samples
    mosdepth_info = mosdepth.get_info(input_dir, low_coverage_threshold=low_coverage_threshold)
    sample_stat_info = stat_info(mosdepth_info, low_coverage_threshold=low_coverage_threshold)

    # Get Variant matrix using for Variant Heatmap
    #
    mutation = []
    variant_matrix_data = []
    samples_variants_info = variants.get_info(input_dir)
    if samples_variants_info:
        df_variants = variants.to_dataframe(samples_variants_info.values())
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
                             depths_data=depths_data,
                             variants_data=variants_data,
                             ref_seq=ref_seq,
                             coverage_stat=sample_stat_info,
                             gene_amplicon_feature_data=gene_amplicon_feature_data,
                             gene_feature_name=gene_feature_name,
                             about_html=about_html,
                             output_html=output_html,
                             region_amplicon_depth_data=region_amplicon_depth_data,
                             variant_matrix_data=variant_matrix_data,
                             mutation=mutation,
                             amplicon=amplicon,
                             gene_feature=gene_feature,
                             low_coverage_threshold=low_coverage_threshold)
