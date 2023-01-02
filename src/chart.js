import {getDatasets} from "./chartOptions/datasets";
import {getXAxes, getYAxes} from "./chartOptions/axes";
import {maxBy, some} from "lodash";
import {getDepthSeries, getRegionAmpliconDepthSeries, getVariantsSeries} from "./chartOptions/series";
import {getGeneFeatureSeries} from "./features";
import {getFluGeneFeature, getSegmentCoords, getYAxisMax} from "./chartOptions/flu/getFluSegmentsInfo";
import {getFluVariantSeries} from "./chartOptions/flu/getFluVariantSeries";
import {getFluPrimerSeries} from "./chartOptions/flu/primerMatches";
import {getTooltips} from "./chartOptions/tooltips";
import {getGrids} from "./chartOptions/grids";
import {getDataZoom} from "./chartOptions/datazoom";


/**
 * Build ECharts options Object from wgscovplot DB object data
 * @param {Object} db - wgscovplot DB object
 * @returns {Object} - ECharts options object
 */
function getCoverageChartOption(db) {
    const {
        echart_features,
        selectedSamples,
        show_amplicons,
        show_genes,
    } = db;
    // TODO: is doubleStrand for finding out if there are gene features being plotted?
    db.doubleStrand = some(echart_features, {value: {strand: -1}});
    const series = [
            ...getDepthSeries(db),
            ...getVariantsSeries(db),
    ];
    if (show_amplicons || show_genes) {
        series.push(...getRegionAmpliconDepthSeries(db));
        series.push(getGeneFeatureSeries({db, index: selectedSamples.length}));
    }
    return {
        title: {},
        dataset: getDatasets(db),
        xAxis: getXAxes(db),
        yAxis: getYAxes(db),
        series: series,
        tooltip: getTooltips(db),
        toolbox: {
            show: "true",
            feature: {
                saveAsImage: {
                    name: "wgscovplot",
                },
                restore: {},
            },
        },
        dataZoom: getDataZoom(selectedSamples.length),
        grid: getGrids(db)
    };
}

/**
 * Get Coverage Chart options
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @returns {Object} - Coverage Chart Option
 */
function getFluCoverageChartOption(db) {
    const {
        selectedSamples,
        segments,
    } = db;
    let chartOptions = {};
    if (selectedSamples.length === 0 || segments.length === 0) {
        return chartOptions; // Plot nothing
    }
    db.segCoords = getSegmentCoords(db);
    let geneFeatureData = getFluGeneFeature(db);
    db.yMax = getYAxisMax(db);
    const totalMaxSegmentLength = maxBy(Object.values(db.segCoords), "end").end;
    let positions = [...Array(totalMaxSegmentLength + 1).keys()];
    positions.shift();
    chartOptions = {
        title: {},
        dataset: getDatasets(db),
        xAxis: getXAxes(db),
        yAxis: getYAxes(db),
        series: [
            ...getDepthSeries(db),
            ...getFluVariantSeries(db),
            ...getFluPrimerSeries(db),
            ...getGeneFeatureSeries(db),
        ],
        tooltip: tooltips(db),
        toolbox: {
            show: "true",
            feature: {
                saveAsImage: {
                    name: "wgscovplot",
                },
                restore: {},
            },
        },
        dataZoom: getDataZoom(selectedSamples.length),
        grid: getGrids(db),
    };
    return chartOptions;

}

export {
    getFluCoverageChartOption,
    getCoverageChartOption,
    getGrids,
    getDataZoom,
    getTooltips,
};