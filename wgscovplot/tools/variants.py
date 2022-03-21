"""VCF and SnpEff/SnpSift parsing functions"""
import logging
import os
import re
from enum import Enum
from operator import itemgetter
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Iterable, Union

import pandas as pd
from pydantic import BaseModel

from wgscovplot.util import try_parse_number, find_file_for_each_sample

logger = logging.getLogger(__name__)

aa_codes = dict(
    ALA='A',
    ARG='R',
    ASN='N',
    ASP='D',
    CYS='C',
    GLU='E',
    GLN='Q',
    GLY='G',
    HIS='H',
    ILE='I',
    LEU='L',
    LYS='K',
    MET='M',
    PHE='F',
    PRO='P',
    SER='S',
    THR='T',
    TRP='W',
    TYR='Y',
    VAL='V',
    TER='*',
)

# list of tuples: [0]: column id; [1]: report column id/name; [2]: column description for report cell comment
variants_cols = [
    ('sample', 'Sample', 'Sample name',),
    ('CHROM', 'Reference Genome', 'Reference genome sequence ID/name'),
    (
        'mutation',
        'Mutation',
        'Mutation found in sample with format '
        '"{reference allele}{reference position}{allele in sample}"'
        ' with predicted amino acid change information in brackets with format '
        '"{gene name}:{reference AA}{gene AA position}{AA change}"'
    ),
    ('POS', 'Position', '1-based nucleotide position in reference sequence'),
    ('REF', 'Reference Allele', 'Nucleotide allele sequence found in reference sequence'),
    ('ALT', 'Alternate Allele', 'Nucleotide allele sequence found in sample'),
    (
        'REF_DP',
        'Reference Allele Depth',
        'Read depth of coverage supporting reference allele at reference position'
    ),
    (
        'ALT_DP',
        'Alternate Allele Depth',
        'Read depth of coverage supporting alternate allele at reference position',
    ),
    ('DP', 'Total Depth', 'Total depth of coverage at reference position',),
    (
        'ALT_FREQ',
        'Alternate Allele Frequency',
        'Observed frequency of alternate allele variant',
    ),
    ('gene', 'Gene', 'Gene name',),
    (
        'impact',
        'Variant Impact',
        'SnpEff estimation of putative impact or deleteriousness of variant '
        '(see https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files)'
    ),
    (
        'effect',
        'Variant Effect',
        'Effect of variant annotated using Sequence Ontology terms, e.g.'
        'for "missense_variant", see http://www.sequenceontology.org/browser/current_release/term/SO:0001583'
        ' where the definition is "A sequence variant, that changes one or more bases, resulting in a '
        'different amino acid sequence but where the length is preserved."'
    ),
    ('aa', 'Amino Acid Change', 'The change in the sample\'s gene amino acid sequence'
                                ' relative to the reference sequence'),
    ('aa_pos', 'Amino Acid Position', 'Position of amino acid change in the reference sequence gene'),
    ('aa_len', 'Gene Amino Acid Length', 'Amino acid length of the reference sequence gene'),
]

variant_summary_cols = [
    (
        'Mutation',
        'Mutation',
        'Mutation found in sample with format '
        '"{reference allele}{reference position}{allele in sample}"'
        ' with predicted amino acid change information in brackets with format '
        '"{gene name}:{reference AA}{gene AA position}{AA change}"'
    ),
    ('n_samples', '# of Samples', 'Number of samples with the mutation.'),
    ('samples', 'Samples', 'List of samples with mutation delimited by semicolon (";")'),
    (
        'min_depth',
        'Min Depth',
        'Minimum depth in all samples that the mutation is observed at.'
    ),
    (
        'max_depth',
        'Max Depth',
        'Maximum depth in all samples that the mutation is observed at.'
    ),
    (
        'mean_depth',
        'Mean Depth',
        'Mean/average depth that the mutation is observed at.'
    ),
    (
        'min_af',
        'Min AF',
        'Minimum alternate allele frequency of mutation in all samples.'
    ),
    (
        'max_af',
        'Max AF',
        'Maximum alternate allele frequency of mutation in all samples.'
    ),
    (
        'mean_af',
        'Mean AF',
        'Mean/average alternate allele frequency of mutation in all samples.'
    ),
    (
        'nt_pos',
        'Nucleotide Position',
        'Nucleotide position of mutation with respect to reference genome.'
    ),
    (
        'aa_pos',
        'Amino Acid Position',
        'Amino acid position of mutation in reference genome gene.'
    )
]

BCFTOOLS_STATS_GLOB_PATTERNS = [
    '**/ivar/**/*AF0.*.bcftools_stats.txt',
    '**/ivar/**/*.bcftools_stats.txt',
    '**/*AF0.*.bcftools_stats.txt',
    '**/*.bcftools_stats.txt',
]

BCFTOOLS_STATS_SAMPLE_NAME_CLEANUP = [
    re.compile(r'\.AF0\.\d+\.bcftools_stats.txt$')
]


class VariantStats(BaseModel):
    sample: str
    n_snp: int
    n_mnp: int
    n_indel: int


class VariantCaller(Enum):
    iVar = 'iVar'
    Longshot = 'Longshot'
    Nanopolish = 'nanopolish'
    Medaka = 'medaka'


VCF_GLOB_PATTERNS = [
    '**/nanopolish/*.pass.vcf.gz',
    '**/ivar/*.vcf.gz',
    '**/*.filt.no_fs.vcf',
    '**/*.longshot.vcf',
    '**/*.vcf',
]

VCF_SAMPLE_NAME_CLEANUP = [
    re.compile(r'(\.pass)?\.vcf(\.gz)?$'),
    re.compile(r'\.AF0\.\d+(\.filt)?'),
    re.compile(r'\.0\.\d+AF(\.filt)?'),
    re.compile(r'\.medaka'),
    re.compile(r'\.longshot'),
    re.compile(r'\.snpeff'),
    re.compile(r'\.no_fs'),
]

SNPSIFT_GLOB_PATTERNS = [
    '**/ivar/**/*.snpSift.table.txt',
    '**/ivar/**/*.snpsift.table.txt',
    '**/ivar/**/*.snpsift.txt',
    '**/*.snpSift.table.txt',
    '**/*.snpsift.table.txt',
    '**/*.snpsift.txt',
]

SNPSIFT_SAMPLE_NAME_CLEANUP = [
    re.compile(r'\.snp[sS]ift\.table\.txt$'),
    re.compile(r'\.snp[sS]ift\.txt$'),
    re.compile(r'\.AF0\.\d+(\.filt)?'),
    re.compile(r'\.0\.\d+AF(\.filt)?'),
]


def vcf_selector(paths: List[Path]) -> Optional[Path]:
    xs = []
    for path in paths:
        variant_caller, df = read_vcf(path)
        xs.append((df.shape[0], path))
    xs.sort(reverse=True)
    try:
        return xs[0][1]
    except KeyError:
        return None


def snpsift_selector(paths: List[Path]) -> Optional[Path]:
    xs = []
    for path in paths:
        df = pd.read_table(path)
        xs.append((df.shape[0], path))
    xs.sort(reverse=True)
    try:
        return xs[0][1]
    except KeyError:
        return None


def read_vcf(vcf_file: Path) -> Tuple[str, pd.DataFrame]:
    """Read VCF file into a DataFrame"""
    gzipped = vcf_file.name.endswith('.gz')
    with os.popen(f'zcat < {vcf_file.absolute()}') if gzipped else open(vcf_file) as fh:
        vcf_cols = []
        variant_caller = ''
        for line in fh:
            if line.startswith('##source='):
                variant_caller = line.strip().replace('##source=', '')
            if line.startswith('##nanopolish'):
                variant_caller = 'nanopolish'
            if line.startswith('##medaka_version'):
                variant_caller = 'medaka'
            if line.startswith('#CHROM'):
                vcf_cols = line[1:].strip().split('\t')
                break
        df = pd.read_table(fh,
                           comment='#',
                           header=None,
                           names=vcf_cols)
        df = df[~df.duplicated(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER'], keep='first')]
    return variant_caller, df


def parse_aa(gene: str,
             ref: str,
             alt: str,
             nt_pos: int,
             aa_pos: int,
             snpeff_aa: str,
             effect: str) -> str:
    if snpeff_aa == '.':
        return f'{ref}{nt_pos}{alt}'
    m = re.match(r'p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)', snpeff_aa)
    if m is None and snpeff_aa.startswith('p.'):
        aa_str = snpeff_aa[2:].upper()
        for aa3, aa1 in aa_codes.items():
            aa_str = aa_str.replace(aa3, aa1)
        if aa_str.endswith('DEL'):
            aa_str = aa_str.replace('DEL', 'del')
        return f'{gene}:{aa_str} ({ref}{nt_pos}{alt})'
    ref_aa, aa_pos_str, alt_aa = m.groups()
    ref_aa = get_aa(ref_aa)

    if effect == 'stop_lost':
        alt_aa = get_aa(alt_aa.replace('ext', ''))
        return f'{gene}:{ref_aa}{aa_pos_str}{alt_aa} ({ref}{nt_pos}{alt} [stop_lost])'
    if effect == 'frameshift_variant':
        return f'{gene}:{ref_aa}{aa_pos_str} ({ref}{nt_pos}{alt} [FRAMESHIFT])'
    if effect == 'conservative_inframe_deletion':
        return f'{gene}:{ref_aa}{aa_pos_str}{alt_aa} ({ref}{nt_pos}{alt})'
    if effect == 'disruptive_inframe_deletion':
        return f'{gene}:{ref_aa}{aa_pos_str}{alt_aa} ({ref}{nt_pos}{alt} [disruptive_inframe_deletion])'

    alt_aa = get_aa(alt_aa)
    if alt_aa == ref_aa:
        return f'{ref}{nt_pos}{alt}'
    return f'{gene}:{ref_aa}{aa_pos_str}{alt_aa} ({ref}{nt_pos}{alt})'


def get_aa(s: str) -> str:
    out = ''
    for i in range(0, len(s), 3):
        aa = s[i: i + 3]
        try:
            aa_code = aa_codes[aa.upper()]
        except KeyError:
            aa_code = aa
        out += aa_code
    return out


def simplify_snpsift(df: pd.DataFrame, sample_name: str) -> Optional[pd.DataFrame]:
    if df.empty:
        return None
    df = df[~df.duplicated(keep='first')]
    field_names = set()
    series = []
    for c in df.columns:
        if c == 'AC':
            df_ac = df[c].str.split(',', n=1, expand=True)
            REF_AC = df_ac[0].astype(int)
            REF_AC.name = 'REF_AC'
            ALT_AC = df_ac[1].astype(int)
            ALT_AC.name = 'ALT_AC'
            AF = ALT_AC / (REF_AC + ALT_AC)
            AF.name = 'AF'
            series += [REF_AC, ALT_AC, AF]
            continue
        if c == 'SR':
            df_ac = df[c].str.split(',', expand=True)
            REF_AC = df_ac[0].astype(int) + df_ac[1].astype(int)
            REF_AC.name = 'REF_AC'
            ALT_AC = df_ac[2].astype(int) + df_ac[3].astype(int)
            ALT_AC.name = 'ALT_AC'
            AF = ALT_AC / (REF_AC + ALT_AC)
            AF.name = 'AF'
            series += [REF_AC, ALT_AC, AF]
            continue
        idx = c.find('[*].')
        if idx > 0:
            new_series_name = c[idx + 4:].lower()
            if new_series_name in field_names:
                continue
            else:
                field_names.add(new_series_name)
            dfc = df[c]
            if dfc.dtype == 'object' and isinstance(dfc.values[0], str):
                new_series = dfc.str.split(',', n=1, expand=True)[0]
            else:
                new_series = dfc.astype('str')
            new_series.name = new_series_name
            series.append(new_series)
        else:
            series.append(df[c])
    df_out = pd.concat(series, axis=1)
    mutation_desc = []
    for row in df_out.itertuples():
        mutation_desc.append(parse_aa(gene=row.gene,
                                      ref=row.REF,
                                      alt=row.ALT,
                                      nt_pos=row.POS,
                                      aa_pos=row.aa_pos,
                                      snpeff_aa=row.aa,
                                      effect=row.effect))

    df_out['mutation'] = mutation_desc
    df_out['sample'] = sample_name
    return df_out


def parse_ivar_vcf(df: pd.DataFrame, sample_name: str = None) -> Optional[pd.DataFrame]:
    if df.empty:
        return None
    if not sample_name:
        sample_name = df.columns[-1] if df.columns[-1] != 'SAMPLE' else None
        if sample_name is None:
            raise ValueError(f'Sample name is not defined for VCF: shape={df.shape}; columns={df.columns}')
    pos_fmt_val = {}
    for row in df.itertuples():
        # An iVar VCF has the DP, total depth, in the INFO field, you cannot sum the alt and ref depths since the ref
        # depth only applies to the first base in the ref allele rather than over the entire allele so the depth can be
        # misleading for deletions
        # See iVar issue: https://github.com/andersen-lab/ivar/issues/86
        infos = parse_vcf_info(row.INFO)
        total_dp = infos['DP']
        ks = row.FORMAT.split(':')
        vs = row[-1].split(':')
        record: Dict[str, Union[float, int, str]] = {k: try_parse_number(v) for k, v in zip(ks, vs)}
        ref_dp = record['REF_DP']
        alt_dp = record['ALT_DP']
        # if the sum of the ref and alt dp does not equal the total dp reported by iVar then recalculate the ref dp
        # since it is likely only reporting the ref dp for the first base of a longer deletion. SNPs should be fine.
        sum_ref_alt_dp = ref_dp + alt_dp
        if sum_ref_alt_dp != total_dp:
            record['REF_DP'] = total_dp - alt_dp
            logger.warning(
                f'iVar VCF for sample "{sample_name}" contains a variant at position {row.POS} where sum of ref allele '
                f'depth ({ref_dp}) and alt allele depth ({alt_dp}) does not equal the total depth ({total_dp}) '
                f'reported by iVar (i.e. {ref_dp} + {alt_dp} = {sum_ref_alt_dp}; {sum_ref_alt_dp} != {total_dp}). Ref '
                f'allele is "{row.REF}", alt allele is "{row.ALT}". iVar alt allele frequency: {record["ALT_FREQ"]}')
        pos_fmt_val[row.POS] = record
        pos_fmt_val[row.POS]['DP'] = total_dp
    df_ivar_info = pd.DataFrame(pos_fmt_val).transpose()
    df_ivar_info.index.name = 'POS'
    df_ivar_info.reset_index(inplace=True)
    df_merge = pd.merge(df, df_ivar_info, on='POS')
    df_merge['sample'] = sample_name
    return df_merge.drop(columns=['ID', 'INFO', 'QUAL', 'FILTER', 'FORMAT', df.columns[-1], 'GT'])


def merge_vcf_snpsift(df_vcf: Optional[pd.DataFrame],
                      df_snpsift: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    if df_snpsift is None and df_vcf is None:
        return None
    if df_vcf is None:
        snpsift_cols = set(df_snpsift.columns)
        return df_snpsift.loc[:, [x for x, _, _ in variants_cols if x in snpsift_cols]]
    if df_snpsift is None:
        vcf_cols = set(df_vcf.columns)
        return df_vcf.loc[:, [x for x, _, _ in variants_cols if x in vcf_cols]]

    df_merge = pd.merge(df_vcf, df_snpsift)
    merged_cols = set(df_merge.columns)
    return df_merge.loc[:, [x for x, _, _ in variants_cols if x in merged_cols]]


def parse_longshot_vcf(df: pd.DataFrame, sample_name: str = None) -> Optional[pd.DataFrame]:
    if df.empty:
        return None
    if not sample_name:
        sample_name = df.columns[-1] if df.columns[-1] != 'SAMPLE' else None
        if sample_name is None:
            raise ValueError(f'Sample name is not defined for VCF: shape={df.shape}; columns={df.columns}')
    pos_info_val = {}
    for row in df.itertuples():
        infos = parse_vcf_info(row.INFO)
        ac_ref, ac_alt = infos['AC']
        infos['REF_DP'] = ac_ref
        infos['ALT_DP'] = ac_alt
        pos_info_val[row.POS] = infos
    df_longshot_info = pd.DataFrame(pos_info_val).transpose()
    df_longshot_info.index.name = 'POS'
    df_longshot_info.reset_index(inplace=True)
    df_merge = pd.merge(df, df_longshot_info, on='POS')
    df_merge['sample'] = sample_name
    df_merge = df_merge[df_merge.DP > 0]
    df_merge['ALT_FREQ'] = df_merge.ALT_DP / df_merge.DP
    cols_to_keep = list({col for col, _, _ in variants_cols} & set(df_merge.columns))
    return df_merge.loc[:, cols_to_keep]


def parse_medaka_vcf(df: pd.DataFrame, sample_name: str = None, min_coverage: int = 10) -> Optional[pd.DataFrame]:
    if df.empty:
        return None
    if not sample_name:
        sample_name = df.columns[-1] if df.columns[-1] != 'SAMPLE' else None
        if sample_name is None:
            raise ValueError(f'Sample name is not defined for VCF: shape={df.shape}; columns={df.columns}')
    pos_info_val = {}
    for row in df.itertuples():
        # if there is no variant INFO available then this isn't the right file
        if row.INFO == '.':
            return None
        infos = parse_vcf_info(row.INFO)
        # no DP INFO? skip this file
        if 'DP' not in infos:
            return None
        if infos['DP'] < min_coverage:
            continue
        ac_ref_fwd, ac_ref_rev, ac_alt_fwd, ac_alt_rev = infos['SR']
        infos['REF_DP'] = ac_ref_fwd + ac_ref_rev
        infos['ALT_DP'] = ac_alt_fwd + ac_alt_rev
        if infos['REF_DP'] + infos['ALT_DP'] < min_coverage:
            continue
        pos_info_val[row.POS] = infos
    if pos_info_val == {}:
        return None
    df_medaka_info = pd.DataFrame(pos_info_val).transpose()
    df_medaka_info.index.name = 'POS'
    df_medaka_info.reset_index(inplace=True)
    df_merge = pd.merge(df, df_medaka_info, on='POS')
    df_merge['sample'] = sample_name
    df_merge = df_merge[df_merge.DP > 0]
    df_merge['ALT_FREQ'] = df_merge.ALT_DP / (df_merge.ALT_DP + df_merge.REF_DP)
    cols_to_keep = list({col for col, _, _ in variants_cols} & set(df_merge.columns))
    return df_merge.loc[:, cols_to_keep]


def parse_nanopolish_vcf(df: pd.DataFrame, sample_name: str = None) -> Optional[pd.DataFrame]:
    if df.empty:
        return None
    if not sample_name:
        sample_name = df.columns[-1] if df.columns[-1] != 'sample' else None
        if sample_name is None:
            raise ValueError(f'Sample name is not defined for VCF: shape={df.shape}; columns={df.columns}')
    pos_info_val = {}
    for row in df.itertuples():
        infos = parse_vcf_info(row.INFO)
        fwd_ref, rev_ref, fwd_alt, rev_alt = infos['StrandSupport']
        allele_count = infos['AlleleCount']
        if allele_count > 1:
            raise NotImplementedError(f'Handling of allele count of {allele_count} is not supported. '
                                      f'Only allele counts of 1 are supported.')
        infos['DP'] = fwd_ref + fwd_alt + rev_ref + rev_alt
        infos['REF_DP'] = fwd_ref + rev_ref
        infos['ALT_DP'] = fwd_alt + rev_alt
        pos_info_val[row.POS] = infos
    df_nanopolish_info = pd.DataFrame(pos_info_val).transpose()
    df_nanopolish_info.index.name = 'POS'
    df_nanopolish_info.reset_index(inplace=True)
    df_merge = pd.merge(df, df_nanopolish_info, on='POS')
    df_merge['sample'] = sample_name
    df_merge = df_merge[df_merge.DP > 0]
    df_merge['ALT_FREQ'] = df_merge.ALT_DP / df_merge.DP
    cols_to_keep = list({col for col, _, _ in variants_cols} & set(df_merge.columns))
    return df_merge.loc[:, cols_to_keep]


def parse_vcf_info(s: str) -> dict:
    out = {}
    for x in s.split(';'):
        if not x:
            continue
        key, val_str = x.split('=', maxsplit=1)
        out[key] = try_parse_number(val_str)
    return out


def get_info(basedir: Path, min_coverage: int = 10) -> Dict[str, pd.DataFrame]:
    sample_vcf = find_file_for_each_sample(basedir=basedir,
                                           glob_patterns=VCF_GLOB_PATTERNS,
                                           sample_name_cleanup=VCF_SAMPLE_NAME_CLEANUP,
                                           single_entry_selector_func=vcf_selector)
    sample_dfvcf = {}
    for sample, vcf_path in sample_vcf.items():
        variant_caller, df_vcf = read_vcf(vcf_path)
        if variant_caller.startswith(VariantCaller.iVar.value):
            df_parsed_ivar_vcf = parse_ivar_vcf(df_vcf, sample)
            if df_parsed_ivar_vcf is not None:
                sample_dfvcf[sample] = df_parsed_ivar_vcf
            else:
                logger.warning(f'Sample "{sample}" has no entries in VCF "{vcf_path}"')
        elif variant_caller.startswith(VariantCaller.Medaka.value):
            df_medaka_vcf = parse_medaka_vcf(df_vcf, sample, min_coverage)
            if df_medaka_vcf is not None:
                sample_dfvcf[sample] = df_medaka_vcf
            else:
                logger.warning(f'Sample "{sample}" has no entries in VCF "{vcf_path}"')
        elif variant_caller.startswith(VariantCaller.Longshot.value):
            df_longshot_vcf = parse_longshot_vcf(df_vcf, sample)
            if df_longshot_vcf is not None:
                sample_dfvcf[sample] = df_longshot_vcf
            else:
                logger.warning(f'Sample "{sample}" has no entries in VCF "{vcf_path}"')
        elif variant_caller.startswith(VariantCaller.Nanopolish.value):
            df_nanopolish_vcf = parse_nanopolish_vcf(df_vcf, sample)
            if df_nanopolish_vcf is not None:
                sample_dfvcf[sample] = df_nanopolish_vcf
            else:
                logger.warning(f'Sample "{sample}" has no entries in VCF "{vcf_path}"')
        else:
            logger.warning(f'Sample "{sample}" VCF file "{vcf_path}" with variant_caller={variant_caller} not supported. Skipping...')

    sample_snpsift = find_file_for_each_sample(basedir=basedir,
                                               glob_patterns=SNPSIFT_GLOB_PATTERNS,
                                               sample_name_cleanup=SNPSIFT_SAMPLE_NAME_CLEANUP,
                                               single_entry_selector_func=snpsift_selector)
    if not sample_snpsift:
        logger.warning(f'No SnpSift tables found in "{basedir}" using glob patterns "{SNPSIFT_GLOB_PATTERNS}"')
    sample_dfsnpsift = {}
    for sample, snpsift_path in sample_snpsift.items():
        df_snpsift = simplify_snpsift(pd.read_table(snpsift_path), sample)
        if df_snpsift is not None:
            sample_dfsnpsift[sample] = df_snpsift
        else:
            logger.warning(f'Sample "{sample}" has no entries in VCF "{snpsift_path}"')
    out = {}
    set_vcf_samples = set(sample_dfvcf.keys())
    set_snpsift_samples = set(sample_dfsnpsift.keys())
    all_samples = set_vcf_samples | set_snpsift_samples
    logger.debug(f'all_samples={len(all_samples)} | '
                 f'vcf only samples={set_vcf_samples - set_snpsift_samples} |'
                 f'snpsift only samples={set_snpsift_samples - set_vcf_samples}')
    for sample in all_samples:
        df_snpsift = sample_dfsnpsift.get(sample, None)
        df_vcf = sample_dfvcf.get(sample, None)
        df_merged = merge_vcf_snpsift(df_vcf, df_snpsift)
        if df_merged is None:
            continue
        out[sample] = df_merged

    return out


def to_dataframe(dfs: Iterable[pd.DataFrame]) -> pd.DataFrame:
    df = pd.concat(list(dfs))
    df.sort_values(['sample', 'POS'], inplace=True)
    df.set_index('sample', inplace=True)
    df.index.name = 'Sample'
    return df.rename(columns={x: y for x, y, _ in variants_cols})


def get_nt_position_int(s: str) -> int:
    """Get the nucleotide position from the mutation string

    >>> get_nt_position_int('A1879G')
    1879
    >>> get_nt_position_int('S:Y145_H146del (TATTACC21992T)')
    21992
    """
    if ':' in s:
        return int(re.sub(r'.*\([AGTC]+(\d+).*', r'\1', s))
    else:
        return int(re.sub(r'[AGTC]+(\d+).*', r'\1', s))


def to_variant_pivot_table(df: pd.DataFrame) -> pd.DataFrame:
    df_vars = df.copy()
    df_vars.reset_index(inplace=True)
    df_pivot = pd.pivot_table(df_vars,
                              index='Sample',
                              columns='Mutation',
                              values='Alternate Allele Frequency',
                              aggfunc='first',
                              fill_value=0.0)
    nt_positions = [get_nt_position_int(x) for x in df_pivot.columns]
    pivot_cols = list(zip(df_pivot.columns,
                          nt_positions))
    pivot_cols.sort(key=itemgetter(1))
    return df_pivot[[x for x, y in pivot_cols]]


def to_summary(df: pd.DataFrame) -> pd.DataFrame:
    df_vars = df.copy()
    df_vars.reset_index(inplace=True)
    logger.debug(f'df_vars columns: {df_vars.columns}')
    df_summary = df_vars.groupby('Mutation', sort=False).agg(
        n_samples=('Sample', 'size'),
        samples=('Sample', lambda x: '; '.join(x)),
        gene=('Gene', 'first'),
        effect=('Variant Effect', 'first'),
        impact=('Variant Impact', 'first'),
        aa=('Amino Acid Change', 'first'),
        min_depth=('Alternate Allele Depth', 'min'),
        max_depth=('Alternate Allele Depth', 'max'),
        mean_depth=('Alternate Allele Depth', lambda x: sum(x) / len(x)),
        min_af=('Alternate Allele Frequency', 'min'),
        max_af=('Alternate Allele Frequency', 'max'),
        mean_af=('Alternate Allele Frequency', lambda x: sum(x) / len(x)),
        nt_pos=('Position', 'first'),
        aa_pos=('Amino Acid Position', 'first')
    )
    df_summary.sort_values('nt_pos', inplace=True)
    return df_summary.rename(columns={x: y for x, y, _ in (variant_summary_cols + variants_cols)})
