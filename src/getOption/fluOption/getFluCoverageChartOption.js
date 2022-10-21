import {getDataSet} from "../getDataSet";
import {getGrids} from "../getGrids";
import {getYAxes} from "../getAxes";
import {getXAxes} from "../getAxes";
import {getDepthSeries} from "../getDepthSeries";
import {getMaxSegmentsLength, getSegmentsInterval} from "./getFluSegmentsInfo";
import {getFluGeneFeature} from "./getFluSegmentsInfo";
import {sum} from "lodash/math";
import {getGeneFeatureSeries} from "../../geneFeatures/getGeneFeatureSeries";
import {getFluTooltips} from "./getFluTooltips";
import {getYAxisMax} from "./getFluSegmentsInfo";
import {getFluVariantSeries} from "./getFluVariantSeries";
import {getFluPrimerSeries} from "./getFluPrimerSeries";
import {getDataZoom} from "../getDataZoom";
import {getToolbox} from "../getToolbox";

/**
 * Get Coverage Chart options
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Object} variants - The object of variants data
 * @param {Object} refSeq - The object of reference sequences of each segment
 * @param {Object} refID - The object of reference ID of each segment
 * @param {Object} lowCoverageRegion - The object of low coverage regions
 * @param {number} lowCoverageThreshold - Low coverage threshold
 * @param {Object} primerData - The object of primer alignment of each segment
 * @param {string} triggerOnType - mousemove or click (tooltips options)
 * @param {boolean} variantSites- whether to show variant sites information (tooltips options)
 * @param {boolean} nonVariantSites - whether to show non-variant sites information (tooltips options)
 * @param {boolean} coverageStatView - whether to show coverage statistics (tooltips options)
 * @param {boolean} infoComparison - whether to compare variants/coverage stat across samples (tooltips options)
 * @param {boolean} showMutation - whether to show mutation below variant sites
 * @param {boolean} showXAxisLabel - whether to show xAxis label
 * @param {boolean} hideOverlapMutation - whether to hide overlapping mutation under variants sites
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 *  variants: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': {}
 *              }
 *          }
 *  refSeq: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': ref_seq
 *              }
 *          }
 *   refID: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': ref_id
 *              }
 *          }
 *   primerData: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': {}
 *              }
 *          }
 */
function getFluCoverageChartOption(samples, segments,
                                   depths, variants,
                                   refSeq, refID,
                                   lowCoverageRegion = {},
                                   lowCoverageThreshold = 10,
                                   primerData = {},
                                   triggerOnType = "mousemove", variantSites = true,
                                   nonVariantSites = false, infoComparison = true,
                                   coverageStatView = false, showMutation = false,
                                   showXAxisLabel = false, hideOverlapMutation = true) {

    let chartOptions = {};
    if (samples.length === 0 || segments.length === 0) {
        return chartOptions; // Plot nothing
    }
    let maxSegmentsLength = getMaxSegmentsLength(samples, segments, depths);
    let segmentsInterval = getSegmentsInterval(maxSegmentsLength);
    let geneFeatureData = getFluGeneFeature(segments, segmentsInterval);
    let yMax = getYAxisMax(samples, segments, depths);
    let positions = [...Array(sum(maxSegmentsLength) + 1).keys()];
    positions.shift();
    chartOptions = {
        title: {},
        dataset: getDataSet(samples, segments, depths, positions),
        xAxis: getXAxes(samples, segments, segmentsInterval, positions.length, true, false, showXAxisLabel ),
        yAxis: getYAxes(samples, "log", yMax, true, false),
        series: [
            ...getDepthSeries(samples, segments, depths, lowCoverageThreshold, segmentsInterval, nonVariantSites),
            ...getFluVariantSeries(samples, segments, depths, variants, segmentsInterval,
                refSeq, variantSites, showMutation, hideOverlapMutation),
            ...getFluPrimerSeries(samples, segments, primerData, segmentsInterval),
            ...getGeneFeatureSeries(geneFeatureData, samples.length, true, false)
        ],
        tooltip: getFluTooltips(samples, segments, depths, variants, refSeq, refID,
            triggerOnType, infoComparison, coverageStatView, lowCoverageThreshold, primerData),
        toolbox: getToolbox(),
        dataZoom: getDataZoom(samples),
        grid: getGrids(samples, true, false, false)
    };
    return chartOptions;

}

export {getFluCoverageChartOption};