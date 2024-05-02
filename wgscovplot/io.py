import logging
from collections import OrderedDict
from io import TextIOWrapper
from pathlib import Path
from typing import Optional, Union

from BCBio import GFF
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from jinja2 import Environment, FileSystemLoader

from wgscovplot.db import DB, SegmentedGenomeDB
from wgscovplot.features import Feature
from wgscovplot.tools import mosdepth
from wgscovplot.util import get_ref_name_bam

logger = logging.getLogger(__name__)

REF_FASTA_GLOB_PATTERN = "*.fa*"
GFF_GLOB_PATTERN = "*.gff"


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


def parse_gff(path: Path) -> list[Feature]:
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
                qs: OrderedDict[str, list[str]] = feature.qualifiers
                feature_name = qs["Name"][0] if "Name" in qs else qs["gbkey"][0]
                out.append(Feature(start=start_pos, end=end_pos, strand=strand, name=feature_name))
    out.sort(key=lambda k: k.start)
    return out


def parse_genbank(gb_handle_or_path: Union[Path, TextIOWrapper]) -> tuple[str, list[Feature]]:
    """Parse sequence and features from Genbank

    Args:
        gb_handle_or_path: Handle or path to Genbank file
    """
    features: list[Feature] = []
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
                qs: OrderedDict[str, list[str]] = seq_feature.qualifiers
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
) -> tuple[Optional[str], Optional[list[Feature]]]:
    try:
        ref_id = mosdepth.get_refseq_id(input_dir)
    except Exception as e:
        logger.error(f"Error getting ref seq ID: {e}")
        ref_id = None
    if not ref_id:
        ref_id = get_refseq_id_from_bam(input_dir)
    if ref_id is None:
        logger.error("Could not find reference sequence ID from BAM file. Exiting.")
        return None, None
    gff_path = find_ref_gff(input_dir, refseq_id=ref_id)
    gene_features = parse_gff(gff_path) if gff_path else None
    fasta_path = find_ref_fasta(input_dir, refseq_id=ref_id)
    ref_seq = read_first_seq_from_fasta(fasta_path) if fasta_path else None
    if ref_seq is None or gene_features is None:
        logger.info(f"Could not find GFF or parse FASTA. Fetching from NCBI for ref seq ID {ref_id}")
        ref_seq, gene_features = fetch_ref_seq_from_ncbi_entrez(ref_id)
    return ref_seq, gene_features


def fetch_ref_seq_from_ncbi_entrez(ref_id: str) -> tuple[Optional[str], Optional[list[Feature]]]:
    if not ref_id:
        logger.info("No ref seq ID provided. Cannot fetch seq from NCBI.")
        return None, None
    logger.info(f'Fetching reference sequence "{ref_id}" from NCBI Entrez')
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ref_id) as handle:
        return parse_genbank(handle)


def read_first_seq_from_fasta(fasta_path: Path) -> str:
    with open(fasta_path) as fh:
        for _, seq in SeqIO.FastaIO.SimpleFastaParser(fh):
            return seq
    return ""


def find_ref_fasta(input_dir: Path, refseq_id: str) -> Optional[Path]:
    for path in input_dir.rglob(REF_FASTA_GLOB_PATTERN):
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line.split()[0]
                    if refseq_id in header:
                        return path
                    else:
                        break
                if line.strip() == "":
                    continue
                else:
                    break
    return None


def find_ref_gff(input_dir: Path, refseq_id: str) -> Optional[Path]:
    gff_paths = list(input_dir.rglob(GFF_GLOB_PATTERN))
    for gff_path in gff_paths:
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                first_col = line.split("\t")[0]
                if refseq_id in first_col:
                    return gff_path
                else:
                    break
    return None


def get_refseq_id_from_bam(input_dir: Path) -> Optional[str]:
    try:
        bam_path = next(input_dir.rglob("*.bam"))
        return get_ref_name_bam(bam_path) if bam_path else None
    except Exception as e:
        logger.error(f"Error getting ref seq ID from BAM file: {e}")
        ref_id = None
    return ref_id
