import logging
from io import TextIOWrapper
from pathlib import Path
from typing import List, Dict, Union, Tuple, OrderedDict, Optional

from BCBio import GFF
from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from jinja2 import Environment, FileSystemLoader
from pydantic import BaseModel

from wgscovplot.db import TemplateDB
from wgscovplot.features import Feature
from wgscovplot.tools import mosdepth

logger = logging.getLogger(__name__)

REF_FASTA_GLOB_PATTERN = '**/genome/snpeff_db/**/*.fa*'
GFF_GLOB_PATTERN = '**/genome/snpeff_db/**/*.gff'


class TemplateHTML(BaseModel):
    about_html: str
    cov_stats_html: str


def write_html_coverage_plot_segment_virus(
        samples_name: List[str],
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
                                        about_html=about_html,
                                        segment_virus=True))


def write_html_coverage_plot(
        db: TemplateDB,
        html: TemplateHTML,
        output_html: Path,
        dev: bool = False,
) -> None:
    render_env = Environment(
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(Path.joinpath(Path(__file__).resolve().parent, "tmpl")),
    )
    template_file = render_env.get_template("wgscovplot.html")
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(
            template_file.render(
                db=db.dict(),
                **html.dict(),
                is_segmented=False,
                dev=dev,
            )
        )


def pydantic_to_dict(x):
    if isinstance(x, list):
        return [pydantic_to_dict(y) for y in x]
    if isinstance(x, dict):
        return {k: pydantic_to_dict(v) for k, v in x.items()}
    if isinstance(x, BaseModel):
        return x.dict()


def parse_gff(path: Path) -> List[Feature]:
    out = []
    with open(path) as fh:
        interest_info = dict(gff_type=["gene", "five_prime_UTR", "three_prime_UTR"])
        rec: SeqRecord
        for rec in GFF.parse(fh, limit_info=interest_info):
            feature: SeqFeature
            for feature in rec.features:
                start_pos = int(feature.location.start) + 1
                end_pos = int(feature.location.end)
                strand = int(feature.location.strand)
                qs: OrderedDict[str, List[str]] = feature.qualifiers
                feature_name = qs['Name'][0] if 'Name' in qs else qs['gbkey'][0]
                out.append(Feature(
                    start_pos=start_pos,
                    end_pos=end_pos,
                    strand=strand,
                    name=feature_name
                ))
    out.sort(key=lambda k: k.start_pos)
    return out


def parse_genbank(
        gb_handle_or_path: Union[Path, TextIOWrapper]
) -> Tuple[str, List[Feature]]:
    """Parse sequence and features from Genbank

    Args:
        gb_handle_or_path: Handle or path to Genbank file
    """
    #features: List[Feature] = []
    features = {}
    skip_feature = {"CDS", "source", "repeat_region", "misc_feature"}
    seq = ''
    for seq_record in SeqIO.parse(gb_handle_or_path, "gb"):
        seq = str(seq_record.seq)
        seq_feature: SeqFeature
        for seq_feature in seq_record.features:
            if seq_feature.type in skip_feature:
                continue
            if seq_feature.type in ["5'UTR", "3'UTR"]:
                feature_name = seq_feature.type
            else:
                qs: OrderedDict[str, List[str]] = seq_feature.qualifiers
                gene = qs.get('gene')
                locus_tag = qs.get('locus_tag')
                if gene:
                    feature_name = gene[0]
                elif locus_tag:
                    feature_name = locus_tag[0]
                else:
                    feature_name = seq_feature.id
            start_pos = int(seq_feature.location.start) + 1
            end_pos = int(seq_feature.location.end)
            strand = int(seq_feature.strand)
            if features.get(strand):
                features[strand].append(Feature(
                    start=start_pos,
                    end=end_pos,
                    strand=strand,
                    name=feature_name
                ))
            else:
                features[strand] = []
                features[strand].append(Feature(
                    start=start_pos,
                    end=end_pos,
                    strand=strand,
                    name=feature_name
                ))
    # sort plus, minus strand coord separately
    if 1 in features.keys():
        features[1].sort(key=lambda k: k.start)
    if -1 in features.keys():
        features[-1].sort(key=lambda k: k.start)
    return seq, features


def get_ref_seq_and_annotation(
        input_dir: Path,
        get_gene_features: bool
) -> Tuple[str, Optional[List[Feature]]]:
    ref_seq = None
    gene_features = None
    if get_gene_features:
        gff_path = find_ref_gff(input_dir)
        gene_features = parse_gff(gff_path) if gff_path else None
    fasta_path = find_ref_fasta(input_dir)
    if fasta_path:
        ref_seq = read_first_seq_from_fasta(fasta_path)
    if ref_seq is None or (gene_features is None and get_gene_features):
        ref_id = mosdepth.get_refseq_id(input_dir)
        ref_seq, gene_features = fetch_ref_seq_from_ncbi_entrez(ref_id)
        if not get_gene_features:
            gene_features = None
    return ref_seq, gene_features


def fetch_ref_seq_from_ncbi_entrez(ref_id: str) -> Tuple[str, List[Feature]]:
    logger.info(f'Fetching reference sequence "{ref_id}" from NCBI Entrez')
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ref_id) as handle:
        return parse_genbank(handle)


def read_first_seq_from_fasta(fasta_path: Path) -> str:
    with open(fasta_path) as fh:
        for _, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
            return seq


def find_ref_fasta(input_dir):
    ref_fasta_path = None
    ref_fasta_paths = list(input_dir.glob(REF_FASTA_GLOB_PATTERN))
    if ref_fasta_paths:
        ref_fasta_path = ref_fasta_paths[0]
    return ref_fasta_path


def find_ref_gff(input_dir: Path) -> Optional[Path]:
    path = None
    gff_paths = list(input_dir.glob(GFF_GLOB_PATTERN))
    if gff_paths:
        path = gff_paths[0]
    return path
