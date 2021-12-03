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
 * @param {Array<Object>} geneFeatureAmpliconData - Array of dictionary geneFeature or amplicon data
 * @param {Array<Object>} ampliconDepthBarData - Array of dictionary geneFeature or amplicon data
 * @param {string} refSeq - Reference seq
 * @param {number} yAxisMax - Max of Y Axis
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {string} geneFeature - whether to plot gene feature or not ("True" or "False")
 * @param {string} amplicon - whether to plot amplicon feature or not ("True" or "False")
 * @param {string} triggerOnType - mousemove or click
 * @param {boolean} isVariantSitesOnly - whether to show tooltips for variant sites only
 * @param {boolean} isVariantComparation - whether to compare variants across samples
 * @returns {Dict[]} - The options for coverage chart
 *
 * The format of data
 * * geneFeatureAmpliconData = [{
 *                    "idx": index,
 *                    "start": start_pos,
 *                    "end": end_pos,
 *                    "level": level,
 *                    "strand": strand,
 *                    "type": "gene_feature or amplicon"}]
 * * ampliconDepthBarData = [{
 *                    "value": [start, end, depth, name]
 *                    "itemStyle": {"color": "skyblue or violet"}}]
 * * samples = ["sample1", "sample2","sample3"]
 * * depths = [
 *       [1, 2, 45, 3, 2, 34, 54, 65, 7, 6, 34, 45, 56, 67, 78, 78],
 *       [1, 2, 45, 0, 0, 9, 15, 65, 7, 6, 20, 8, 4, 15, 100, 102],
 *       [12, 12, 425, 3, 2, 10, 12, 9, 7, 6, 1, 45, 45, 67, 87, 97]
 *   ]
 * * variants = [
 *       [{sample: 'sample1', CHROM: 'MN908947.3',mutation: 'C14408T(orf1ab:L4715L)',POS: 14408},
 *        {sample: 'sample1', CHROM: 'MN908947.3',mutation: 'C14408T(orf1ab:L4715L)',POS: 29300}],
 *       [{sample: 'sample2', CHROM: 'MN908947.3',mutation: 'C15480A(orf1ab:P5072H)',POS: 10000},
 *        {sample: 'sample2', CHROM: 'MN908947.3',mutation: 'C15480A(orf1ab:P5072H)',POS: 25300}]
 *    ]
 */
function getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData,refSeq,
                                yAxisMax, samples, depths, variants,
                                geneFeature= "False", amplicon= "False",
                                triggerOnType= "mousemove", isVariantSitesOnly = false,
                                isVariantComparation = false) {

    var positions = [...Array(refSeq.length + 1).keys()];
    positions.shift();
    var chartOptions = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, positions.length, geneFeature, amplicon),
        yAxis: getYAxes(samples, "log", yAxisMax, geneFeature, amplicon),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths, refSeq),
            ...getAmpliconDepthSeries(samples, ampliconDepthBarData, amplicon),
            ...getGeneFeatureSeries(geneFeatureAmpliconData, samples.length, geneFeature, amplicon)
        ],
        tooltip: getTooltips(samples, depths, variants, refSeq, triggerOnType, isVariantSitesOnly, isVariantComparation),
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