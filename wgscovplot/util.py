import contextlib
import logging
import re
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Optional, Union

import numpy as np
import pandas as pd

ListOfStrOrPattern = Union[list[str], list[re.Pattern[str]], list[Union[str, re.Pattern[str]]]]

logger = logging.getLogger(__name__)

NT_MAP = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


def find_file_for_each_sample(
    basedir: Path,
    glob_patterns: list[str],
    sample_name_cleanup: Optional[ListOfStrOrPattern] = None,
    single_entry_selector_func: Optional[Callable] = None,
) -> Mapping[str, Path]:
    sample_files = defaultdict(list)
    for glob_pattern in glob_patterns:
        for p in basedir.glob(glob_pattern):
            sample = extract_sample_name(p.name, remove=sample_name_cleanup)

            sample_files[sample].append(p)
    sample_file = {}
    for sample, files in sample_files.items():
        if single_entry_selector_func:
            sample_file[sample] = single_entry_selector_func(files)
        else:
            # select first file if no selector func specified
            sample_file[sample] = files[0]
    return sample_file


def select_most_recent_file(files: list[Path]) -> Path:
    return sorted(files, key=lambda x: x.stat().st_mtime, reverse=True)[0]


def extract_sample_name(
    filename: str,
    remove: Optional[ListOfStrOrPattern] = None,
) -> str:
    if not remove:
        remove = [
            ".pass",
            ".mapped",
            ".trim",
            ".ivar_trim",
            ".mkD",
            ".sorted",
            ".bam",
            ".flagstat",
            ".stats",
            ".txt",
            ".idxstats",
            ".depths.tsv",
            "-depths.tsv",
            ".tsv",
            ".mosdepth",
            ".per-base",
            ".bed",
            ".gz",
        ]
    out = filename
    for x in remove:
        if isinstance(x, str):
            out = out.replace(x, "")
        elif isinstance(x, re.Pattern):
            out = x.sub("", out)
        else:
            logger.warning(
                f'Not sure how to use "{x}" (type={type(x)}) for removing '
                f'non-sample name text from "{filename}". Skipping "{x}". Output="{out}".'
            )
    return out


def get_col_widths(
    df: pd.DataFrame, index: bool = False, offset: int = 2, max_width: Optional[int] = None, include_header: bool = True
) -> Iterable[int]:
    """Calculate column widths based on column headers and contents"""
    if index:
        idx_max = max([len(str(s)) for s in df.index.values] + [len(str(df.index.name))]) + offset
        if max_width:
            idx_max = min(idx_max, max_width)
        yield idx_max
    for c in df.columns:
        # get max length of column contents and length of column header
        max_width_cells = df[c].astype(str).str.len().max() + 1
        if include_header:
            width = np.max([max_width_cells, len(c) + 1]) + offset
        elif isinstance(c, str):
            col_words = c.split()
            col_words.sort(key=len, reverse=True)
            max_word_size = int(len(col_words[-1]) * 1.25 + 1)
            width = np.max([max_width_cells, max_word_size]) + offset
        else:
            width = max_width_cells + offset
        if max_width:
            width = min(width, max_width)
        yield width


def get_row_heights(df: pd.DataFrame, idx: int, offset: int = 0, multiplier: int = 15):
    """Calculate row heights"""
    # get max number of newlines in the row
    newline_count = np.max(df.loc[idx, :].astype(str).str.count("\n").max())
    newline_count = max(newline_count, 1)
    height = newline_count * multiplier + offset
    logger.debug(f'idx="{idx}" height={height} newline_count={newline_count}')
    return height


def try_parse_number(s: str) -> Any:
    if "," in s:
        xs = s.split(",")
        return [try_parse_number(x) for x in xs]
    with contextlib.suppress(ValueError):
        return int(s)
    with contextlib.suppress(ValueError):
        return float(s)
    return s


def expand_degenerate_bases(seq: str) -> Iterable[str]:
    for x in product(*(NT_MAP[nt] for nt in seq)):
        yield "".join(x)


def overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 < start2 < end1 or start1 < end2 < end1


def get_ref_name_bam(path: Path) -> str:
    import pysam

    bam = pysam.AlignmentFile(path)
    logger.info(f"BAM: {bam}")
    ref_name = bam.get_reference_name(0)
    logger.info(f"{ref_name=}")
    if not ref_name:
        for ref_name in bam.references:
            if ref_name:
                return ref_name
    return ref_name
