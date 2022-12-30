import {getDatasets} from "./chartOptions/datasets";
import {getXAxes, getYAxes} from "./chartOptions/axes";
import {toTableHtml} from "./util";
import {find} from "lodash/collection";
import {genomeCoverage, meanCoverage, medianCoverage} from "./stats";
import {isNil, maxBy, some} from "lodash";
import {getDepthSeries, getRegionAmpliconDepthSeries, getVariantsSeries} from "./chartOptions/series";
import {getGeneFeatureSeries} from "./features";
import {getFluGeneFeature, getSegmentCoords, getYAxisMax} from "./chartOptions/flu/getFluSegmentsInfo";
import {getFluVariantSeries} from "./chartOptions/flu/getFluVariantSeries";
import {tooltips} from "./chartOptions/flu/tooltips";
import {getCoverageStatComparison, getVariantComparison} from "./chartOptions/tooltips";
import {getFluPrimerSeries} from "./chartOptions/flu/primerMatches";
import {toFloat32Array} from "./util";

/**
 * Get ECharts Grid options objects given selected samples (and segments) and if gene/amplicon/segment features are being shown
 * @param {WgsCovPlotDB} db
 */
function getGrids(db) {
    const {
        selectedSamples,
        show_genes,
        show_amplicons,
        doubleStrand,
        padTop = 4.0,
        heightOffset = 6.0,
    } = db;
    let featureHeight;
    // TODO: make magic numbers function optional params
    if (show_amplicons && show_genes) {
        featureHeight = 15.0;
    } else if (show_amplicons || show_genes) {
        featureHeight = 6.0;
    } else {
        featureHeight = -5.0; // make subplot full
    }
    // TODO: find out what doubleStrand is
    featureHeight = (doubleStrand && featureHeight > 0) ? (featureHeight + 6.0) : featureHeight;
    let subPlotHeight = 90.0 - featureHeight;
    let grids = [];
    let sampleHeight = (subPlotHeight - padTop) / selectedSamples.length - heightOffset;
    for (let idx = 0; idx < selectedSamples.length; idx++) {
        grids.push({
            show: true,
            height: sampleHeight.toFixed(1) + "%",
            top: ((sampleHeight + heightOffset) * idx + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    if (show_amplicons || show_genes) {
        grids.push({
            show: false,
            height: featureHeight.toFixed(1) + "%",
            top: ((sampleHeight + heightOffset) * selectedSamples.length + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    return grids;
}

/**
 * Get ECharts DataZoom options
 * @param {number} numSamples - Number of selected samples
 * @returns {[Object, Object]}
 */
function getDataZoom(numSamples) {
    let xAxisIndex = [...Array(numSamples + 1).keys()];
    return [
        {
            type: "inside",
            filterMode: "none",
            xAxisIndex: xAxisIndex,
            zoomLock: false,
        },
        {
            show: true,
            filterMode: "none",
            xAxisIndex: xAxisIndex,
            type: "slider",
            zoomLock: false,
            showDataShadow: false
        },
    ];
}

function tooltipFormatter({db}) {
    return function (params) {
        const {
            selectedSamples,
            depths,
            chart,
            crossSampleComparisonInTooltips,
            variants,
            ref_seq,
            low_coverage_threshold,
            showCovStatsInTooltips,
        } = db;
        let output = "";
        let [{
            axisIndex,
            axisValue: position
        }] = params;
        if (axisIndex >= selectedSamples.length) {
            return output;
        }
        let sample = selectedSamples[axisIndex];
        let sampleDepths = toFloat32Array(depths[sample]);
        let depth = sampleDepths[position - 1];
        let [dataZoom] = chart.getOption().dataZoom;
        let zoomStart = Math.floor(dataZoom.startValue);
        let zoomEnd = Math.floor(dataZoom.endValue);
        let positionRows = [];
        let coverageStatRows = [];
        const isVariantBar = params.find(element => {
            return element.componentSubType === "bar";

        });
        if (isVariantBar) {
            if (crossSampleComparisonInTooltips) {
                positionRows = getVariantComparison({
                    db,
                    position,
                    sampleInFocus: sample,
                });
            } else {
                positionRows = [
                    ["Position", position.toLocaleString()],
                    ["Coverage Depth", depth.toLocaleString()],
                ];
                let foundObj = find(variants[sample], {"POS": position.toString()}, 0);
                if (!isNil(foundObj)) {
                    for (const [key, value] of Object.entries(foundObj)) {
                        if (key !== "POS" && key !== "sample") {
                            positionRows.push([key, value]);
                        }
                    }
                }
            }
        } else {
            positionRows = [
                ["Pos", position.toLocaleString()],
                ["Coverage Depth", depth.toLocaleString()],
            ];
            positionRows.push(["Sequence", ref_seq[position - 1]]);
        }
        if (positionRows.length) {
            output += `<h5>Sample: ${sample}</h5>`;
            output += toTableHtml({
                headers: ["Position Info", ""],
                rows: positionRows,
            });
            if (showCovStatsInTooltips) {
                if (crossSampleComparisonInTooltips) {
                    coverageStatRows = getCoverageStatComparison({
                        db,
                        start: zoomStart,
                        end: zoomEnd,
                        position
                    });
                } else {
                    let meanCov = meanCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
                    let medianCov = medianCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
                    let genomeCov = genomeCoverage(sampleDepths, zoomStart, zoomEnd, low_coverage_threshold).toFixed(2);
                    coverageStatRows = [
                        [
                            "Range",
                            `${zoomStart.toLocaleString()} - ${zoomEnd.toLocaleString()}`,
                        ],
                        ["Mean Coverage", `${meanCov}X`],
                        ["Median Coverage", `${medianCov}X`],
                        [`Genome Coverage (>= ${low_coverage_threshold}X)`, `${genomeCov}%`],
                    ];
                }
                output += toTableHtml(
                    {
                        headers: ["Coverage View Stats", ""],
                        rows: coverageStatRows,
                    });
            }
        }
        return output;
    };
}

/**
 * Define options for tooltips
 * @param {Object} db - wgscovplot DB object
 * @returns {Array<Object>}
 */
function getTooltips(db) {
    const {
        tooltipTriggerOn = "mousemove",
    } = db;
    return [
        {
            trigger: "axis",
            enterable: true,
            triggerOn: tooltipTriggerOn,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            position: "cursor",
            axisPointer: {
                type: "line"
            },
            formatter: tooltipFormatter({db}),
        },
    ];
}

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