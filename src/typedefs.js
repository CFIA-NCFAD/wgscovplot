/**
 * @typedef {Object<string, Object<string, number[]>>} SampleSegmentDepths
 */

/**
 * @typedef {Object<string, number[]>} SampleDepths
 */

/**
 * @typedef {Object<string, number[]>} SampleDepths
 */


/**
 * @typedef {Object<string, string>} VariantCall
 */

/**
 * @typedef {Object<string, VariantCall[]>} SampleVariantCalls
 */

/**
 * @typedef {Object<string, Object<string, VariantCall[]>>} SampleSegmentVariantCalls
 */

/**
 * @typedef {Object} AmpliconDepth
 * @property {{color: string}} itemStyle - style of bar chart elements
 * @property {[number, number, number, string]} value - array containing amplicon start, end, depth and name
 */

/**
 * @typedef {Object<string, AmpliconDepth[]>} SampleAmpliconDepths
 */

/**
 * @typedef {Object<string, {maxLength: number, start: number, end: number}>} SegmentCoords
 */

/**
 * @typedef {Object} MosdepthInfo
* @property {number} genome_coverage
* @property {string} low_coverage_coords
* @property {number} low_coverage_threshold
* @property {number} max_depth
* @property {number} mean_coverage
* @property {number} median_coverage
* @property {number} n_low_coverage
* @property {number} n_zero_coverage
* @property {number} ref_seq_length
* @property {string} sample
* @property {string} zero_coverage_coords
 */

/**
 * @typedef {Object<string, MosdepthInfo>} SampleMosdepthInfo
 */

/**
 * @typedef {Object<string, Object<string, MosdepthInfo>>} SampleSegmentMosdepthInfo
 */

/**
 * @typedef WgsCovPlotDB
 * @type {object}
 * @property {string | null} ref_seq
 * @property {SampleDepths | SampleSegmentDepths} depths
 * @property {SampleVariantCalls | SampleSegmentVariantCalls} variants
 * @property {SampleAmpliconDepths | null} amplicon_depths
 * @property {SampleMosdepthInfo | SampleSegmentMosdepthInfo} mosdepth_info
 * @property {Object[] | null} echart_features
 * @property {Object} chart - ECharts object for coverage plots
 * @property {Object} variantHeatmap - ECharts object for variant calling heatmap
 * @property {number[]} positions - array of positions from 1 to length of reference sequence
 * @property {string[]} samples - array of all samples
 * @property {string[]} selectedSamples - Selected sample names
 * @property {boolean} doubleStrand - ??? gene/amplicon features on both plus and minus strands?
 * @property {number} low_coverage_threshold - Low coverage threshold
 * @property {boolean} crossSampleComparisonInTooltips - show cross sample comparison in tooltips?
 * @property {boolean} fixedTooltipPosition - tooltip appears in fixed position?
 * @property {boolean} hideOverlappingVariantLabels - hide overlapping variant labels
 * @property {boolean} show_amplicons - Show amplicon features?
 * @property {boolean} show_genes - Show gene features?
 * @property {boolean} showCovStatsInTooltips - show coverage stats in tooltips
 * @property {boolean} showGeneLabels - show gene feature labels
 * @property {boolean} showNonVariantSiteTooltips - default = false
 * @property {boolean} showVariantLabels - Show variant site labels below each cov subplot x-axis?
 * @property {boolean} showVariantSiteTooltips - Focus tooltips on variant sites?
 * @property {boolean} showDataZoomSlider - Show ECharts DataZoom slider?
 * @property {boolean} tooltipEnabled - Show tooltips for coverage plots?
 * @property {number} heightOffset - Cov subplot height offset as % (default 6.0)
 * @property {number} padTop - Padding to top coordinate of each cov subplot as % (default 4.0)
 * @property {boolean} showLowCovRegions - Highlight low coverage regions in coverage plots?
 * @property {number} showLowCovRegionsOpacity - Opacity of areas highlighting low coverage regions (default 0.5)
 * @property {number} yMax - y-axis max value. Can be adjusted by user, but also dependent on samples being shown.
 * @property {string} scaleType - Cov plot y-axis scaling, "linear" or "log"
 * @property {string} leftMargin - cov plot left margin (default "4.0%")
 * @property {string} rightMargin - cov plot right margin (default "4.0%")
 * @property {string} tooltipTriggerOn - what event to trigger show tooltip on (default "mousemove")
 * @property {string[] | null} segments - list of virus segment names if segmented virus being shown
 * @property {string[] | null} selectedSegments - list of selected virus segment names
 * @property {SegmentCoords | null} segCoords - plotting coordinates for virus segments
 */

/**
 * @typedef Elements
 * @type {object}
 * @property {HTMLDivElement} $chart
 * @property {HTMLDivElement} $variantHeatmap
 * @property {HTMLInputElement} $chartHeightInput
 * @property {HTMLOutputElement} $chartHeightOutput
 * @property {HTMLInputElement} $chartTopInput
 * @property {HTMLOutputElement} $chartTopOutput
 * @property {HTMLInputElement} $chartLeftInput
 * @property {HTMLOutputElement} $chartLeftOutput
 * @property {HTMLInputElement} $chartRightInput
 * @property {HTMLOutputElement} $chartRightOutput
 * @property {HTMLInputElement} $genefeatureHeightInput
 * @property {HTMLOutputElement} $genefeatureHeightOutput
 * @property {HTMLSelectElement} $renderEnv
 * @property {Object} $selectedSamples
 * @property {Object} $selectedGeneFeatures
 * @property {HTMLInputElement} $toggleDarkMode
 * @property {HTMLInputElement} $toggleGeneLabel
 * @property {HTMLInputElement} $toggleAmplicons
 * @property {HTMLInputElement} $toggleShowVariantLabels
 * @property {HTMLInputElement} $startPos
 * @property {HTMLInputElement} $endPos
 * @property {HTMLButtonElement} $btnShowRegion
 * @property {HTMLButtonElement} $btnResetView
 * @property {HTMLInputElement} $ymax
 * @property {HTMLButtonElement} $btnSetYMax
 * @property {HTMLSelectElement} $selectYScale
 * @property {HTMLInputElement} $lowThreshold
 * @property {HTMLButtonElement} $btnSetLowCovThreshold
 * @property {HTMLInputElement} $toggleTooltip
 * @property {HTMLInputElement} $toggleShowAmplicon
 * @property {HTMLInputElement} $toggleHideOverlappingVariantLabels
 * @property {HTMLInputElement} $toggleXAxisLabel
 * @property {HTMLInputElement} $toggleShowVariantSiteTooltips
 * @property {HTMLInputElement} $toggleShowNonVariantSiteTooltips
 * @property {HTMLInputElement} $toggleShowLowCovRegions
 * @property {HTMLInputElement} $toggleFixedTooltipPosition
 * @property {HTMLButtonElement} $btnShowSelectedFeatures
 * @property {HTMLInputElement} $toggleTooltipCovStats
 * @property {HTMLInputElement} $toggleShowTooltipOnClick
 * @property {HTMLInputElement} $toggleTooltipVariantComparison
 * @property {HTMLInputElement} $toggleDataZoomSlider
 * @property {HTMLButtonElement} $btnResetMargins
 * @property {HTMLButtonElement} $btnShowZoomRegion
 * @property {HTMLButtonElement} $btnResetView
 */