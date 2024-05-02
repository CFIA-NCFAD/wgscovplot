import base64
import logging
from collections import defaultdict
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from pydantic import BaseModel

from wgscovplot.colors import AmpliconColour
from wgscovplot.features import Feature
from wgscovplot.util import find_file_for_each_sample, get_ref_name_bam

logger = logging.getLogger(__name__)

SAMPLE_NAME_CLEANUP = [
    ".genome.per-base.bed.gz",
    ".amplicon.per-base.bed.gz",
    ".per-base.bed.gz",
    ".genome.regions.bed.gz",
    ".amplicon.regions.bed.gz",
    ".regions.bed.gz",
    "-depths.tsv",
    ".trim",
    ".ivar_trim",
    ".sorted",
    ".mkD",
    ".topsegments.csv",
    ".bam",
]

PER_BASE_PATTERNS = [
    "**/mosdepth/**/*.genome.per-base.bed.gz",
    "**/mosdepth/**/*.per-base.bed.gz",
    "**/mosdepth/**/*-depths.tsv",
    # fallback to parsing BAM files
    "**/*.trim*.bam",
    "**/*.bam",
]

TOP_REFERENCE_PATTERNS = [
    "**/reference_sequences/**/*.topsegments.csv",
]

REGIONS_PATTERNS = [
    "**/mosdepth/**/*.amplicon.regions.bed.gz",
    "**/mosdepth/**/*.regions.bed.gz",
]


class MosdepthDepthInfo(BaseModel):
    sample: str
    n_zero_coverage: int
    zero_coverage_coords: str
    low_coverage_threshold: int = 5
    n_low_coverage: int
    low_coverage_coords: str
    genome_coverage: float
    mean_coverage: float
    median_coverage: float
    ref_seq_length: int
    max_depth: int


def get_samples_name(basedir: Path, is_genome_segmented: bool) -> list:
    glob_patterns = TOP_REFERENCE_PATTERNS if is_genome_segmented else PER_BASE_PATTERNS
    sample_beds = find_file_for_each_sample(
        basedir, glob_patterns=glob_patterns, sample_name_cleanup=SAMPLE_NAME_CLEANUP
    )
    out = list(sample_beds.keys())
    return sorted(out)


def read_mosdepth_bed(p: Path) -> pd.DataFrame:
    return pd.read_table(p, header=None, names=["genome", "start_idx", "end_idx", "depth"])


def read_mosdepth_region_bed(p: Path) -> pd.DataFrame:
    return pd.read_table(p, header=None, names=["genome", "start_idx", "end_idx", "amplicon", "depth"])


def get_interval_coords(depths: np.ndarray, threshold: int = 0) -> str:
    """Get coordinates of intervals where depth is zero or below threshold if specified.

    >>> get_interval_coords(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    '1-10'
    >>> get_interval_coords(np.array([0,1,2,3,4,0,0,0,1,2,3,4,5,6,0,0]))
    '1; 6-8; 15-16'
    >>> get_interval_coords(np.array([0,1,2,3,4,0,0,0,1,2,3,4,5,6,0,0]), threshold=2)
    '1-2; 6-9; 15-16'

    Args:
        depths: Coverage depths array
        threshold: Minimum depth threshold

    Returns:
        Coordinates of intervals where depth is zero or below threshold
    """
    mask = depths == 0 if threshold == 0 else depths < threshold
    below = np.where(mask)[0]
    coords: list[list[int]] = []
    for x in below:
        if coords:
            last = coords[-1][-1]
            if x == last + 1:
                coords[-1].append(x)
            else:
                coords.append([x])
        else:
            coords.append([x])
    return "; ".join([f"{xs[0] + 1}-{xs[-1] + 1}" if xs[0] != xs[-1] else f"{xs[0] + 1}" for xs in coords])


def get_genome_coverage(depths: np.ndarray, low_coverage_threshold: int = 5) -> float:
    """Calculate genome coverage as a fraction of positions with depth >= low_coverage_threshold

    >>> get_genome_coverage(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    0.0
    >>> get_genome_coverage(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 5]))
    0.1
    >>> get_genome_coverage(np.array([10, 9, 8, 7, 6, 5, 4, 3, 2, 1]))
    0.6

    Args:
        depths: Coverage depths array
        low_coverage_threshold: Minimum depth threshold

    Returns:
        Genome coverage as a fraction of positions with depth >= low_coverage_threshold
    """
    return 1.0 - (np.count_nonzero(depths < low_coverage_threshold) / len(depths))


def depth_array(df: pd.DataFrame) -> np.ndarray:
    arr = np.zeros(df.end_idx.max(), dtype=np.uint16)
    for row in df.itertuples():
        arr[row.start_idx : row.end_idx] = row.depth
    return arr


def get_refseq_id(basedir: Path) -> str:
    sample_paths = find_file_for_each_sample(
        basedir, glob_patterns=PER_BASE_PATTERNS, sample_name_cleanup=SAMPLE_NAME_CLEANUP
    )
    refseq_id = ""
    for path in sample_paths.values():
        if path.suffix == ".bam":
            logger.info(f"Trying to get ref seq id from BAM at {path}")
            return get_ref_name_bam(path)
        df = read_mosdepth_bed(path)
        refseq_id = df["genome"][0] if df is not None and not df.empty else ""
        if not refseq_id:
            return refseq_id
    return refseq_id


def get_amplicon_depths(basedir: Path) -> dict[str, list]:
    sample_path = find_file_for_each_sample(
        basedir, glob_patterns=REGIONS_PATTERNS, sample_name_cleanup=SAMPLE_NAME_CLEANUP
    )
    out = defaultdict(list)
    sample = None
    try:
        for sample, path in sample_path.items():
            if path.suffix == ".bam":
                break
            df = read_mosdepth_region_bed(path)
            for row in df.itertuples():
                pool_id = int(row.amplicon.split("_")[-1])
                color = AmpliconColour.pool2 if pool_id % 2 == 0 else AmpliconColour.pool1
                out[sample].append(
                    {
                        "value": [row.start_idx, row.end_idx, row.depth, row.amplicon],
                        "itemStyle": {"color": color},
                    }
                )
        return dict(out)
    except Exception as e:
        logger.error(e, exc_info=True)
        logger.warning(f"{sample} No Region Amplicon Depth Found")
        return {}


def get_region_amplicon(basedir: Path) -> list[Feature]:
    sample_path = find_file_for_each_sample(
        basedir, glob_patterns=REGIONS_PATTERNS, sample_name_cleanup=SAMPLE_NAME_CLEANUP
    )
    try:
        for path in sample_path.values():
            if path.suffix == ".bam":
                break
            df_amplicon = pd.read_table(path, names=["reference", "start", "end", "amplicon", "depth"], header=None)
            return [
                Feature(
                    name=row.amplicon,
                    start=row.start,
                    end=row.end,
                )
                for row in df_amplicon.itertuples()
            ]
    except Warning:
        logger.warning("No Region Amplicon Found")
    return []


def get_info(
    basedir: Path, low_coverage_threshold: int = 5
) -> tuple[dict[str, MosdepthDepthInfo], dict[str, np.ndarray]]:
    sample_beds = find_file_for_each_sample(
        basedir, glob_patterns=PER_BASE_PATTERNS, sample_name_cleanup=SAMPLE_NAME_CLEANUP
    )
    out = {}
    sample_depths = {}
    for sample, path in sample_beds.items():
        logger.info(f"{sample=} {path=}")
        if path.suffix == ".bam":
            arr = parse_bam_depths(path)
        else:
            df = read_mosdepth_bed(path)
            arr = depth_array(df)
        sample_depths[sample] = arr
        mean_cov = arr.mean()
        median_cov = pd.Series(arr).median()
        depth_info = MosdepthDepthInfo(
            sample=sample,
            low_coverage_threshold=low_coverage_threshold,
            n_low_coverage=(arr < low_coverage_threshold).sum(),
            n_zero_coverage=(arr == 0).sum(),
            zero_coverage_coords=get_interval_coords(arr, threshold=0),
            low_coverage_coords=get_interval_coords(arr, threshold=low_coverage_threshold),
            genome_coverage=get_genome_coverage(arr, low_coverage_threshold),
            mean_coverage=mean_cov,
            median_coverage=median_cov,
            ref_seq_length=len(arr),
            max_depth=arr.max(initial=0),
        )
        out[sample] = depth_info
    return out, sample_depths


def get_base64_encoded_depth_arrays(
    sample_depths: Union[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]
) -> Union[dict[str, str], dict[str, dict[str, str]]]:
    """Encode depth arrays as base64 strings

    Instead of dumping a list of numbers to a JSON list, the float32 array will be base64 encoded so
    that it can be decoded into a Float32Array in JS. The base64 encoding will compress the numbers
    significantly (~60% of the size of dumping list to JSON).

    ```python
    import base64
    import numpy as np
    import json

    t = np.arange(30000, dtype=np.float32)
    print(len(base64.b64encode(t)))
    # 160000

    print(json.dumps(t.tolist()))
    # 258890
    ```

    The array can be converted to a JS Float32Array with:

    ```js
    // base64 encoded 32-bit float array
    b64 = "AAAAAAAAgD8AAABAAABAQAAAgEAAAKBAAADAQAAA4..."
    f32a = new Float32Array(new Uint8Array([...atob(b64)].map(c => c.charCodeAt(0))).buffer)
    // to regular JS array
    arr = Array.from(f32a)
    ```

    Args:
        sample_depths: Dict of sample name (to segment name) to depths array

    Returns:
        Dict of sample name (to segment name) to 32-bit float depth arrays encoded with base64
    """
    first_value = next(iter(sample_depths.values()))
    if isinstance(first_value, dict):
        return {
            sample: {segment: to_base64(arr) for segment, arr in segment_depths.items()}
            for sample, segment_depths in sample_depths.items()
        }
    return {sample: to_base64(arr) for sample, arr in sample_depths.items()}


def to_base64(x: np.ndarray) -> str:
    return base64.b64encode(x).decode("utf-8")


def parse_bam_depths(bam_path: Path) -> np.ndarray:
    """Parse BAM file for depths of coverage across first reference genome using pysam

    Args:
        bam_path: Path to BAM file

    Returns:
        Coverage depths
    """
    import pysam

    bam = pysam.AlignmentFile(bam_path)
    reference_name = bam.references[0]
    reference_length = bam.get_reference_length(reference_name)
    depths = np.zeros(reference_length, dtype=np.int32)
    for pileupcolumn in bam.pileup():
        depth = sum(1 for pileupread in pileupcolumn.pileups if not (pileupread.is_del or pileupread.is_refskip))
        depths[pileupcolumn.reference_pos] = depth
    return depths
