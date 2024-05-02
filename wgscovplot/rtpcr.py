from collections import defaultdict
from pathlib import Path

import edlib
from Bio import SeqIO
from pydantic import BaseModel

from wgscovplot.util import expand_degenerate_bases


class PrimerMatch(BaseModel):
    name: str
    cigar: str
    edit_distance: int
    start: int
    end: int
    other_locations: str
    query_aligned: str
    matched_aligned: str
    target_aligned: str


def flu_rtpcr_matches(
    primer_seq_path: Path, sample_seg_seq: dict[str, dict[str, str]], edit_distance_threshold: int
) -> defaultdict[str, defaultdict[str, list[PrimerMatch]]]:
    """Find real-time PCR primer/probe matches to each segment of each sample using Edlib."""
    out: defaultdict[str, defaultdict[str, list[PrimerMatch]]] = defaultdict(lambda: defaultdict(list))
    primer_seq_record = [(record.id, record.seq) for record in SeqIO.parse(open(primer_seq_path), "fasta")]
    for sample, seg_seq in sample_seg_seq.items():
        for seg, seq in seg_seq.items():
            if len(seq):
                for seq_id, primer_seq in primer_seq_record:
                    alns = []
                    for pseq in expand_degenerate_bases(primer_seq):
                        aln = edlib.align(query=pseq, target=seq, mode="HW", task="path", k=edit_distance_threshold)
                        if aln["editDistance"] != -1:  # keep sequence has alignment result
                            aln["query_seq"] = pseq
                            alns.append(aln)
                    if alns:
                        best_aln = min(alns, key=lambda x: x["editDistance"])
                        other_locations = []
                        aln_loc = best_aln["locations"]
                        start, end = aln_loc[0]
                        for loc_start, loc_end in aln_loc[1:]:
                            # convert to 1-based index
                            other_locations.append((loc_start + 1, loc_end + 1))
                        nice_aln = edlib.getNiceAlignment(best_aln, best_aln["query_seq"], seq)
                        match = PrimerMatch(
                            name=seq_id,
                            cigar=best_aln["cigar"],
                            edit_distance=best_aln["editDistance"],
                            start=start,
                            end=end,
                            other_locations=", ".join(map(str, other_locations)),
                            query_aligned=nice_aln["query_aligned"],
                            matched_aligned=nice_aln["matched_aligned"],
                            target_aligned=nice_aln["target_aligned"],
                        )
                        out[sample][seg].append(match)
    return out
