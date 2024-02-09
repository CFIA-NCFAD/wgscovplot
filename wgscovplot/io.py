import logging
from io import TextIOWrapper
from pathlib import Path
from typing import List, Optional, OrderedDict, Tuple, Union

from BCBio import GFF
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from jinja2 import Environment, FileSystemLoader

from wgscovplot.db import DB, SegmentedGenomeDB
from wgscovplot.features import Feature
from wgscovplot.tools import mosdepth

logger = logging.getLogger(__name__)

REF_FASTA_GLOB_PATTERN = "**/genome/**/*.fa*"
GFF_GLOB_PATTERN = "**/genome/**/*.gff"


def write_html_coverage_plot(
    db: Union[DB, SegmentedGenomeDB],
    output_html: Path,
) -> None:
    """Write HTML coverage plot"""
    root = Path(__file__).resolve().parent.parent
    web_dir = root / "web"
    # autoescape=True, # maybe?
    render_env = Environment(  # noqa: S701
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(web_dir),
    )
    template_file = render_env.get_template("wgscovplot.html.j2")
    with open(output_html, "w+", encoding="utf-8") as fout:
        fout.write(
            template_file.render(
                db=db.model_dump(),
            )
        )


def parse_gff(path: Path) -> List[Feature]:
    out = []
    with open(path) as fh:
        interest_info = {"gff_type": ["gene", "five_prime_UTR", "three_prime_UTR"]}
        rec: SeqRecord
        for rec in GFF.parse(fh, limit_info=interest_info):
            feature: SeqFeature
            for feature in rec.features:
                start_pos = int(feature.location.start) + 1
                end_pos = int(feature.location.end)
                strand = int(feature.location.strand)
                qs: OrderedDict[str, List[str]] = feature.qualifiers
                feature_name = qs["Name"][0] if "Name" in qs else qs["gbkey"][0]
                out.append(Feature(start=start_pos, end=end_pos, strand=strand, name=feature_name))
    out.sort(key=lambda k: k.start)
    return out


def parse_genbank(gb_handle_or_path: Union[Path, TextIOWrapper]) -> Tuple[str, List[Feature]]:
    """Parse sequence and features from Genbank

    Args:
        gb_handle_or_path: Handle or path to Genbank file
    """
    features: List[Feature] = []
    skip_feature = {"CDS", "source", "repeat_region", "misc_feature"}
    seq = ""
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
                gene = qs.get("gene")
                locus_tag = qs.get("locus_tag")
                if gene:
                    feature_name = gene[0]
                elif locus_tag:
                    feature_name = locus_tag[0]
                else:
                    feature_name = seq_feature.id
            start_pos = int(seq_feature.location.start) + 1
            end_pos = int(seq_feature.location.end)
            strand = int(seq_feature.strand)
            features.append(Feature(start=start_pos, end=end_pos, strand=strand, name=feature_name))
    features.sort(key=lambda f: f.start)
    return seq, features


def get_ref_seq_and_annotation(
    input_dir: Path,
) -> Tuple[str, Optional[List[Feature]]]:
    gff_path = find_ref_gff(input_dir)
    gene_features = parse_gff(gff_path) if gff_path else None
    fasta_path = find_ref_fasta(input_dir)
    ref_seq = read_first_seq_from_fasta(fasta_path) if fasta_path else None
    if ref_seq is None or gene_features is None:
        ref_id = mosdepth.get_refseq_id(input_dir)
        ref_seq, gene_features = fetch_ref_seq_from_ncbi_entrez(ref_id)
    return ref_seq, gene_features


def fetch_ref_seq_from_ncbi_entrez(ref_id: str) -> Tuple[str, List[Feature]]:
    logger.info(f'Fetching reference sequence "{ref_id}" from NCBI Entrez')
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ref_id) as handle:
        return parse_genbank(handle)


def read_first_seq_from_fasta(fasta_path: Path) -> str:
    with open(fasta_path) as fh:
        for _, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
            return seq
    return ""


def find_ref_fasta(input_dir: Path) -> Optional[Path]:
    ref_fasta_paths = list(input_dir.glob(REF_FASTA_GLOB_PATTERN))
    return ref_fasta_paths[0] if ref_fasta_paths else None


def find_ref_gff(input_dir: Path) -> Optional[Path]:
    gff_paths = list(input_dir.glob(GFF_GLOB_PATTERN))
    return gff_paths[0] if gff_paths else None
