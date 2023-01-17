from itertools import chain
from pathlib import Path
from typing import List, Set, Dict

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pydantic.main import BaseModel

from wgscovplot.tools.mosdepth import TOP_REFERENCE_PATTERNS, SAMPLE_NAME_CLEANUP
from wgscovplot.util import find_file_for_each_sample


class SampleSegmentRef(BaseModel):
    sample: str
    segment: str
    ref_seq: str
    ref_id: str

'''
def find_ref_seqs(
        basedir: Path,
        ref_ids: Set[str]
) -> Dict[str, str]:
    """Try to find sequences in FASTA files

    Args:
        basedir: directory to search for FASTA files
        ref_ids: sequence IDs to look for

    Returns:
        Dict of ref seq ID to nucleotide sequence
    """
    refs = {}
    for fasta_path in chain(basedir.rglob('*.fa'), basedir.rglob('*.fasta')):
        with open(fasta_path) as f:
            for header, seq in SimpleFastaParser(f):
                seq_id, _ = header.split(' ', maxsplit=1)
                if seq_id in ref_ids:
                    refs[seq_id] = seq
    return refs
'''


def read_top_references_table(basedir: Path) -> pd.DataFrame:
    return pd.read_csv(
        basedir,
        header=0,
        names=[
            'sample',
            'segment',
            'ref_id',
            'blastn_bitscore',
            'ref_sequence_id',
        ]
    )


def get_sample_top_references(basedir: Path) -> Dict[str, pd.DataFrame]:
    """Find files sample.topsegments.csv for each sample

    File is located in reference_sequences/sample/sample.topsegments.csv
    """
    out = {}
    sample_top_references = find_file_for_each_sample(
        basedir,
        glob_patterns=TOP_REFERENCE_PATTERNS,
        sample_name_cleanup=SAMPLE_NAME_CLEANUP,
    )
    for sample, top_references_path in sample_top_references.items():
        out[sample] = read_top_references_table(top_references_path)
    return out


'''
def get_sample_top_references(basedir: Path) -> List[SampleSegmentRef]:
    """Find files sample.topsegments.csv for each sample

    File is located in reference_sequences/sample/sample.topsegments.csv
    """
    sample_file = find_file_for_each_sample(
        basedir,
        glob_patterns=TOP_REFERENCE_PATTERNS,
        sample_name_cleanup=SAMPLE_NAME_CLEANUP,
    )
    print (sample_file)
    # TODO: concat all top ref tables if more than one
    # the top ref IDs could also be gotten from BAM files which should always be present
    # this would be more reliable than looking at Mosdepth files, but those should always
    # be output
    df = pd.concat([read_top_references_table(p) for _, p in sample_file.items()])
    ref_ids = set(df['ref_id'])
    refs = find_ref_seqs(basedir, ref_ids)
    return [
        SampleSegmentRef(
            sample=r.sample,
            segment=r.segment,
            ref_id=r.ref_id,
            ref_seq=refs.get(r.ref_id, None),
        ) for r in df.itertuples()]
'''
