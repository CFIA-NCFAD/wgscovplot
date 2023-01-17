import {getDatasets} from "./chartOptions/datasets";
import {getXAxes, getYAxes} from "./chartOptions/axes";
import {isNil, maxBy, some, sum, values} from "lodash";
import {getDepthSeries, getRegionAmpliconDepthSeries, getVariantsSeries} from "./chartOptions/series";
import {getGeneFeatureSeries} from "./features";
import {getFluGeneFeature, getSegmentCoords} from "./chartOptions/flu/segmentInfo";
import {getSegmentVariantSeries} from "./chartOptions/flu/segmentVariantSeries";
import {getSegmentPrimerSeries} from "./chartOptions/flu/primerMatches";
import {getTooltips} from "./chartOptions/tooltips";
import {getGrids} from "./chartOptions/grids";
import {getDataZoom} from "./chartOptions/datazoom";
import {getMaxSegmentLength} from "./chartOptions/flu/segmentInfo";
import {getSegmentTooltips} from "./chartOptions/flu/segmentTooltips";

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
    let dataset
    let series;
    if (db.segment_virus === true){
        let maxSegmentLength = getMaxSegmentLength(db);
        db.segCoords = getSegmentCoords(db, maxSegmentLength);
        db.echart_features = getFluGeneFeature(db);
        db.positions = [...Array(sum(values(maxSegmentLength)) + 1).keys()];
        db.positions.shift();
        dataset = getDatasets(db);
        series = [
            ...getDepthSeries(db),
            ...getSegmentVariantSeries(db),
            ...getSegmentPrimerSeries(db)
        ];
        series.push(getGeneFeatureSeries({db, index: selectedSamples.length}));
    } else {
        // TODO: is doubleStrand for finding out if there are gene features being plotted?
        db.doubleStrand = some(echart_features, {value: {strand: -1}});
        dataset = getDatasets(db);
        series = [
            ...getDepthSeries(db),
            ...getVariantsSeries(db),
        ];
        if (show_amplicons || show_genes) {
            series.push(...getRegionAmpliconDepthSeries(db));
            series.push(getGeneFeatureSeries({db, index: selectedSamples.length}));
        }
    }
    return {
        title: {},
        dataset: dataset,
        xAxis: getXAxes(db),
        yAxis: getYAxes(db),
        series: series,
        tooltip: db.segment_virus === false ? getTooltips(db) : getSegmentTooltips(db),
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

export {
    getCoverageChartOption,
    getGrids,
    getDataZoom,
    getTooltips,
};