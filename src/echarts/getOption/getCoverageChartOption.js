import {getDatasets} from "./getDataSets";
import {getXAxes, getYAxes} from "./getAxes";
import {getDepthSeries} from "./getDepthSeries";
import {getVariantsSeries} from "./getVariantSeries";
import {getAmpliconDepthSeries} from "./getAmpliconDepthSeries";
import {getGeneFeatureSeries} from "../geneFeatures/getGeneFeaturesSeries";
import {getGrids} from "./getGrids";
import {getTooltips} from "./getTooltips";


/**
 * Define all options for coverage chart
 * @param {Array<Dict[]>} geneFeatureAmpliconData - Array of dictionary geneFeature or amplicon data
 * @param {Array<Dict[]>} ampliconDepthBarData - Array of dictionary geneFeature or amplicon data
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @param {string} geneFeature - whether to plot gene feature or not ("True" or "False")
 * @param {string} amplicon - whether to plot amplicon feature or not ("True" or "False")
 * @returns {Dict[]}
 *  * geneFeatureAmpliconData = [{
 *                    "idx": index,
 *                    "start": start_pos,
 *                    "end": end_pos,
 *                    "level": level,
 *"                   "strand": strand,
 *                    "type": "gene_feature or amplicon"}]
 *  * ampliconDepthBarData = [{
 *                    "value": [start, end, depth, name]
 *                    "itemStyle": {"color": "skyblue or violte"}}]
 *
 */
function getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData,
                                samples, depths, variants,
                                geneFeature, amplicon) {
    var chartOptions = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, refSeqLength, geneFeature, amplicon),
        yAxis: getYAxes(samples, "log", maxDepth, geneFeature, amplicon),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths),
            ...getAmpliconDepthSeries(samples, ampliconDepthBarData),
            ...getGeneFeatureSeries(geneFeatureAmpliconData, samples.length, geneFeature, amplicon)
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
                xAxisIndex: [...Array(samples.length + 1).keys()],
                zoomLock: false,
            },
            {
                show: true,
                filterMode: "none",
                xAxisIndex: [...Array(samples.length + 1).keys()],
                type: "slider",
                zoomLock: false,
            },
        ],
        grid: getGrids(samples, geneFeature, amplicon)
    };
    return chartOptions;
}

export {getCoverageChartOption};