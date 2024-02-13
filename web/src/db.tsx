import type {ECharts, EChartsType} from "echarts";

export interface SampleSegmentDepths {
  [key: string]: {
    [key: string]: number[];
  }
}

export interface SampleDepths {
  [key: string]: number[];
}

export interface ECColorArg {
  data: {
    0: number;
  }
}

export interface PrimerMatches {
  [key: string]: {
    [key: string]: PrimerMatchesInfo[]
  }
}

export interface VariantCall {
  [key: string]: string | number;

  POS: string | number;
}

export interface SampleVariantCalls {
  [key: string]: VariantCall[];
}

export interface SampleSegmentVariantCalls {
  [key: string]: {
    [key: string]: VariantCall[];
  }
}

export interface AmpliconDepth {
  /** style of bar chart elements */
  itemStyle: {
    color: string;
  }
  /** array containing amplicon start, end, depth and name */
  value: [number, number, number, string];
}

export interface SampleAmpliconDepths {
  [key: string]: AmpliconDepth[];
}

export interface SegmentCoords {
  [key: string]: {
    maxLength: number;
    start: number;
    end: number;
  }
}

export interface MaxSegmentLength {
    [key: string]: number
}

export interface SegmentRef {
  [key: string]: {
    [key: string] : string
  }
}

export interface MosdepthInfo {
  genome_coverage: number;
  low_coverage_coords: string;
  low_coverage_threshold: number;
  max_depth: number;
  mean_coverage: number;
  median_coverage: number;
  n_low_coverage: number;
  n_zero_coverage: number;
  ref_seq_length: number;
  sample: string;
  zero_coverage_coords: string;
}

export interface PrimerMatchesInfo {
  cigar: string,
  edit_distance: number,
  end: number,
  matched_aligned: string,
  name: string,
  other_locations: string
  query_aligned: string,
  start: number,
  target_aligned: string
}

export interface SampleMosdepthInfo {
  [key: string]: MosdepthInfo;
}

export interface SampleSegmentMosdepthInfo {
  [key: string]: {
    [key: string]: MosdepthInfo;
  }
}

export interface ntColor {
  [key: string]: string;
  A: string;
  C: string;
  G: string;
  T: string;
}

export interface ChartOptions {
  showVariants: boolean;
  variantBarWidth: number;
  showFeatures: boolean;
  geneLabelTextSize: number;
  geneLabelDistance: number;
  geneLabelRotation: number;
  variantLabelsRotation: number;
  coordsLabelsRotation: number;
  lowCovThresholdLineWidth: number;
  lowCovThresholdLineColour: string;
  lowCovColour: string;
  /** Selected sample names */
  selectedSamples: string[];
  /** list of selected virus segment names */
  selectedSegments?: string[];
  /** show coverage stats in tooltips */
  showCovStatsInTooltips: boolean;
  /** Show ECharts DataZoom slider? */
  showDataZoomSlider: boolean;
  /** show gene feature labels */
  showGeneLabels: boolean;
  /** Highlight low coverage regions in coverage plots? */
  showLowCovRegions: boolean;
  /** Opacity of areas highlighting low coverage regions (default 0.5) */
  showLowCovRegionsOpacity: number;
  /** Show variant site labels below each cov subplot x-axis? */
  showVariantLabels: boolean;
  showXAxisLabel: boolean;
  showLowCoverageCoords: boolean;
  /** y-axis max value. Can be adjusted by user, but also dependent on samples being shown. */
  yMax: number;
  covColour: string;
  sidebarCollapsed: boolean;

  /** Show tooltips for coverage plots? */
  tooltipEnabled: boolean;
  /** what event to trigger show tooltip on (default "mousemove") */
  tooltipTriggerOn: string;

  /** show cross sample comparison in tooltips? */
  crossSampleComparisonInTooltips: boolean;

  /** tooltip appears in fixed position? */
  fixedTooltipPosition: boolean;
  /** Cov subplot height offset as % (default 6.0) */
  heightOffset: number;
  /** hide overlapping variant labels */
  hideOverlappingVariantLabels: boolean;
  /** cov plot left margin (default "4.0%") */
  leftMargin: number;
  lowCoverageOpacity: number;
  /** Low coverage threshold */
  low_coverage_threshold: number;
  /** cov plot right margin (default "4.0%") */
  rightMargin: number;
  /** Padding to top coordinate of each cov subplot as % (default 4.0) */
  padTop: number;
  /** Cov plot y-axis scaling, "linear" or "log" */
  scaleType: "value" | "log";
  featurePlotHeightScaling: number;
  subplotTitleFontSize: number;
  subplotTitleColour: string;
  startPos: number;
  endPos: number;
  renderer: "canvas" | "svg";
  darkMode: boolean;
  ntColor: ntColor;
}

export interface Table {
  /** array of column names */
  headers: string[];
  /** array of arrays of row data */
  rows: (string | number)[][];
}

export interface TooltipOptions {
  showTooltip: boolean;
  variantSitesOnly: boolean;
  dragging: boolean;
  /** show tooltips for coverage plots? */
  show: boolean;
  top: number;
  left: number;
  x: number;
  y: number;
  sample: string;
  position: number;
  depth: number;
  tables: Table[];
}

export interface ECFeatureValue {
  start: number;
  end: number;
  idx: number;
  level: number;
  rotate: number;
  type: "gene" | "amplicon" | "segment";
  strand: number
}

export interface ECFeature {
  itemStyle: {
    color: string;
  }
  name: string;
  value: ECFeatureValue;
}

export interface WgsCovPlotDB {
  activePage: string;
  /** Sample amplicon depths if amplicon sequencing was used to generate the sequencing data. */
  amplicon_depths?: SampleAmpliconDepths;
  /** Sample coverage depths by position or sample and segment coverage depths by position for segmented viruses */
  depths: SampleDepths | SampleSegmentDepths;
  /** Are there genetic features on both plus and minus strands? */
  doubleStrand: boolean;
  /** ECharts features */
  echart_features: ECFeature[];
  /** mosdepth/coverage depth info */
  mosdepth_info: SampleMosdepthInfo | SampleSegmentMosdepthInfo;

  /** Array of position numbers from 1 to length of reference sequences or
   * the sum of max genome segment lengths for each selected genome segment.
   */
  positions: number[];
  /** Reference nucleotide sequence for non-segmented viruses. */
  ref_seq?: string;
  /** Reference nucleotide sequences for segmented viruses organized by sample name and segment name. */
  segments_ref_seq?: SegmentRef;
  segments_ref_id?: SegmentRef;

  /** Real-time PCR (rtPCR) primer and probe match data from Edlib analysis if
   * rtPCR primers and probes were specified.
   */
  primer_matches?: PrimerMatches

  /** All sample names, not to be confused with selectedSamples */
  samples: string[];

  /** plotting coordinates for virus segments */
  segCoords?: SegmentCoords;
  /** list of virus segment names if segmented virus being shown */
  segments?: string[];

  /** Max length among samples of 1 segment**/
  maxSegmentLength?: MaxSegmentLength

  /** Show amplicon features? */
  show_amplicons?: boolean;
  /** Show primer matches? */
  show_primer_matches?: boolean;
  /** Show gene features? */
  show_genes: boolean;
  /** ECharts object for variant calling heatmap */
  variantHeatmap?: ECharts | EChartsType;
  /** variant calls */
  variants: SampleVariantCalls | SampleSegmentVariantCalls;

  // eslint-disable-next-line
  chart: any;
  // eslint-disable-next-line
  heatMapChart: any;
  chartOptions: ChartOptions;
  tooltipOptions: TooltipOptions;
  // eslint-disable-next-line
  summaryInfo: any
}

// eslint-disable-next-line
export const defaultDB: WgsCovPlotDB = {
  echart_features: [], heatMapChart: undefined, summaryInfo: undefined,
  activePage: "chart",
  chart: null,
  depths: {},
  doubleStrand: false,
  mosdepth_info: {},
  positions: [],
  ref_seq: "",
  samples: [],
  show_amplicons: false,
  show_primer_matches: false,
  show_genes: true,
  variantHeatmap: undefined,
  variants: {},
  tooltipOptions: {
    showTooltip: true,
    variantSitesOnly: true,
    dragging: false,
    show: false,
    top: 100,
    left: 100,
    x: 0,
    y: 0,
    sample: "",
    position: 0,
    depth: 0,
    tables: [],
  },
  chartOptions: {
    showVariants: true,
    variantBarWidth: 2,
    showFeatures: true,
    geneLabelTextSize: 14,
    geneLabelDistance: 10,
    geneLabelRotation: 0,
    variantLabelsRotation: -45,
    coordsLabelsRotation: 30,
    lowCovThresholdLineWidth: 1,
    lowCovThresholdLineColour: "#910000",
    sidebarCollapsed: false,
    covColour: "#6a6a6a",
    crossSampleComparisonInTooltips: false,
    fixedTooltipPosition: false,
    heightOffset: 6.0,
    hideOverlappingVariantLabels: true,
    leftMargin: 3.0,
    lowCoverageOpacity: 0.3,
    low_coverage_threshold: 10,
    padTop: 3,
    rightMargin: 2.0,
    scaleType: "log",
    tooltipEnabled: false,
    tooltipTriggerOn: "mousemove",
    selectedSamples: [],
    selectedSegments: [],
    showCovStatsInTooltips: true,
    showDataZoomSlider: true,
    showGeneLabels: true,
    showLowCovRegions: false,
    showLowCovRegionsOpacity: 0.5,
    showVariantLabels: false,
    showXAxisLabel: false,
    showLowCoverageCoords: false,
    yMax: 1000,
    lowCovColour: "#ffff00",
    featurePlotHeightScaling: 100,
    subplotTitleFontSize: 16,
    subplotTitleColour: "#232323",
    startPos: 0,
    endPos: 0,
    renderer: "canvas",
    darkMode: false,
    ntColor: {
      A: "#ea5e48",
      C: "#eaca48",
      G: "#6ad82b",
      T: "#2b87d8",
    }
  }
}

export interface ECFormatterFeature {
  componentType: string;
  componentSubType: string;
  componentIndex: number;
  seriesType: string;
  seriesIndex: number;
  seriesId: string;
  seriesName: string;
  name: string;
  dataIndex: number;
  data: ECFeature;
  value: ECFeatureValue;
  color: string;
  marker: string;
}