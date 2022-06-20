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

/**
 * Get Coverage Chart options
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of segments name
 * @param {Object} depths - Object of depths array
 * @param {String} scaleType - The scale set for yAXis either 'log' or 'linear'
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getFluCoverageChartOption(samples, segments, depths, scaleType="log") {

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
        xAxis: getFluXAxes(samples, segments, position.length, true, segmentsRange),
        yAxis: getYAxes(samples, scaleType, yMax, true, false),
        series: [
            ...getFluDepthSeries(samples, segments),
            ...getGeneFeatureSeries(geneFeatureData, samples.length, true, false)
        ],
        tooltip: getFluTooltips(samples, segments, depths, segmentsRange),
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