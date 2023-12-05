from collections import defaultdict
from pathlib import Path
from typing import Dict, DefaultDict

import edlib
from Bio import SeqIO

from wgscovplot.util import expand_degenerate_bases


def flu_rtpcr_matches(
        primer_seq_path: Path,
        sample_seg_seq: Dict[str, Dict[str, str]],
        edit_distance_threshold: int
) -> DefaultDict[str, DefaultDict[str, Dict]]:
    """Find real-time PCR primer/probe matches to each segment of each sample using Edlib."""
    out = defaultdict(lambda: defaultdict(list))
    primer_seq_record = [(record.id, record.seq) for record in SeqIO.parse(open(primer_seq_path), 'fasta')]
    for sample, seg_seq in sample_seg_seq.items():
        for seg, seq in seg_seq.items():
            if len(seq):
                for idx, (seq_id, primer_seq) in enumerate(primer_seq_record):
                    expanded_primer_seqs = expand_degenerate_bases(primer_seq)
                    alns = []
                    for pseq in expanded_primer_seqs:
                        aln = edlib.align(query=pseq,
                                          target=seq,
                                          mode="HW",
                                          task="path",
                                          k=edit_distance_threshold)
                        aln['query_seq'] = pseq
                        if aln['editDistance'] != -1:  # keep sequence has alignment result
                            alns.append(aln)
                    if alns:
                        best_aln = min(alns, key=lambda x: x['editDistance'])
                        other_locations = []
                        aln_loc = best_aln['locations']
                        start, end = aln_loc[0]
                        for loc_start, loc_end in aln_loc[1:]:
                            # convert to 1-based index
                            if loc_start is not None:
                                loc_start = loc_start + 1
                            if loc_end is not None:
                                loc_end = loc_end + 1
                            other_locations.append((loc_start, loc_end))
                        nice_aln = edlib.getNiceAlignment(best_aln, best_aln['query_seq'], seq)
                        match = dict(
                            name=seq_id,
                            cigar=best_aln['cigar'],
                            edit_distance=best_aln['editDistance'],
                            start=start,
                            end=end,
                            other_locations=", ".join(map(str, other_locations))
                        )
                        out[sample][seg].append({**match, **nice_aln})
    return out