import {getDatasets} from "./getDataSets";
import {getXAxes, getYAxes} from "./getAxes";
import {getDepthSeries} from "./getDepthSeries";
import {getVariantsSeries} from "./getVariantSeries";
import {getAmpliconDepthSeries} from "./getAmpliconDepthSeries";
import {getGeneFeatureSeries} from "../geneFeatures/getGeneFeaturesSeries";
import {getGrids} from "./getGrids";
import {getTooltips} from "./getTooltips";

/**
 * A closure is used keep track the current number samples in the chart and is used for toggling slider
 * @type {function(): number}
 */
var gridLength = (function (){
    var len;
    return function (){
        return len;
    }
})();

/**
 * Define all options for coverage chart
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @returns {Dict[]}
 */
function getCoverageChartOption(samples, depths, variants) {
    gridLength = samples.length;
    var chartOptions = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, refSeqLength),
        yAxis: getYAxes(samples, "log", maxDepth),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths),
            ...getAmpliconDepthSeries(samples),
            ...getGeneFeatureSeries(gridLength)
        ],
        tooltip: getTooltips(samples, depths, variants),
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
                xAxisIndex: [...Array(gridLength + 1).keys()],
                zoomLock: false,
            },
            {
                show: true,
                filterMode: "none",
                xAxisIndex: [...Array(gridLength + 1).keys()],
                type: "slider",
                zoomLock: false,
            },
        ],
        grid: getGrids(samples)
    };
    return chartOptions;
}

export {getCoverageChartOption};