import {getFluDatasets} from "./getFluDataSets";
import {getGrids} from "../getGrids";
import {getYAxes} from "../getAxes";
import {getFluXAxes} from "./getFluAxes";
import {getFluDepthSeries} from "./getFluDepthSeries";
import {getMaxSegmentsLength, getSegmentsRange} from "./getFluSegmentsInfo";
import {getFluGeneFeature} from "./getFluSegmentsInfo";
import {sum} from "lodash/math";
import {getGeneFeatureSeries} from "../../geneFeatures/getGeneFeaturesSeries";
import {getFluTooltips} from "./getFluTooltips";
import {getYAxisMax} from "./getFluSegmentsInfo";
import {getFluVariantSeries} from "./getFluVariantSeries";
import {getFluMarkAreaSeries} from "./getFluMarkAreaSeries";

/**
 * Get Coverage Chart options
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Object} variants - The object of variants data
 * @param {Object} refSeq - The object of reference sequences of each segment
 * @param {Object} refID - The object of reference ID of each segment
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
 */
function getFluCoverageChartOption(samples, segments,
                                   depths, variants,
                                   refSeq, refID,
                                   triggerOnType="mousemove", variantSites=true,
                                   nonVariantSites=false, infoComparison=true,
                                   coverageStatView=false, showMutation=false,
                                   showXAxisLabel=false, hideOverlapMutation=true)
{

    let chartOptions = {};
    if (samples.length == 0 || segments.length == 0){
        return chartOptions; // Plot nothing
    }
    let maxSegmentsLength = getMaxSegmentsLength(samples, segments, depths);
    let segmentsRange = getSegmentsRange(maxSegmentsLength);
    let geneFeatureData = getFluGeneFeature(segments, segmentsRange);
    let yMax = getYAxisMax(samples, segments, depths);
    let position = [...Array(sum(maxSegmentsLength) + 1).keys()];
    position.shift();
    chartOptions = {
        title: {},
        dataset: getFluDatasets(samples, segments, depths, position),
        xAxis: getFluXAxes(samples, segments, position.length, showXAxisLabel, segmentsRange),
        yAxis: getYAxes(samples, "log", yMax, true, false),
        series: [
            ...getFluDepthSeries(samples, segments, nonVariantSites),
            ...getFluVariantSeries(samples, segments, depths, variants, segmentsRange,
                refSeq, variantSites, showMutation, hideOverlapMutation),
            ...getGeneFeatureSeries(geneFeatureData, samples.length, true, false)
        ],
        tooltip: getFluTooltips(samples, segments, depths, variants, refSeq, refID,
                                triggerOnType, infoComparison, coverageStatView),
        toolbox: {
            show: "true",
            feature: {
                dataView: {
                    readOnly: false,
                },
                restore: {},
                saveAsImage: {
                    name: "wgscovplot",
                },
            },
        },
        dataZoom: [
            {
                type: "inside",
                filterMode: "none",
                xAxisIndex: [...Array(samples.length + 1).keys()],
                zoomLock: false,
            },
            {
                show: true,
                filterMode: "none",
                xAxisIndex: [...Array(samples.length + 1).keys()],
                type: "slider",
                zoomLock: false,
                showDataShadow: false
            },
        ],
        grid: getGrids(samples, true, false, false)
    };
    return chartOptions;

}

export {getFluCoverageChartOption};