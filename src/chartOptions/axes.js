import {FEATURE_PLOT_PROPS} from "../util";
import {isNil, maxBy} from "lodash";

/**
 * Custom xAxis label (for segment virus)
 * @param {number} value - xAxis value
 * @param {string[]} segments - An array of segments names
 * @param {SegmentCoords} segCoords
 * @returns {string}
 */
function getCustomXAxisLabel(value, segments, segCoords) {
    for (let segment of segments) {
        let {start, end} = segCoords[segment];
        if (value >= start && value <= end) {
            let pos = value - start + 1
            return `${segment} - ${pos.toLocaleString()}`;
        }
    }
    return ""
}

/**
 * Define options for x Axis
 * @param {WgsCovPlotDB} db
 * @returns {Array<Object>}
 */
function getXAxes(db) {
    const {
        selectedSamples,
        selectedSegments = null,
        segCoords = null,
        positions,
        show_genes,
        show_amplicons,
        showXAxisLabel = true,
    } = db;
    let formatter = {};
    if (!isNil(selectedSegments) && selectedSegments.length > 0) {
        formatter.formatter = function (value) {
            return getCustomXAxisLabel(value, selectedSegments, segCoords)
        }
    }
    let axes = [];
    for (let i = 0; i < selectedSamples.length; i++) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: positions.length,
            axisLabel: {
                show: showXAxisLabel,
                interval: "auto",
                ...formatter,
            }
        });
    }
    if (show_amplicons || show_genes) {
        axes.push({
            type: "value",
            gridIndex: selectedSamples.length,
            min: 1,
            max: positions.length,
            axisLabel: {
                interval: "auto",
                ...formatter,
            },
        });
    }
    return axes;
}

/**
 * Define options for Y axis
 * @param {WgsCovPlotDB} db
 * @returns {Object[]}
 */
function getYAxes(db) {
    const {
        selectedSamples,
        scaleType,
        yMax,
        show_genes,
        show_amplicons
    } = db;
    let axes = [];
    for (let [i, sample] of selectedSamples.entries()) {
        axes.push({
            type: scaleType,
            gridIndex: i,
            name: sample,
            nameTextStyle: {
                fontStyle: "normal",
                fontWeight: "bolder",
                fontSize: 12
            },
            nameLocation: "end",
            nameRotate: 0.01,
            min: scaleType === "log" ? 1 : 0,
            max: yMax,
            minorSplitLine: {
                show: true,
            },
        });
    }
    if (show_amplicons || show_genes || db.segment_virus) {
        axes.push({
            max: FEATURE_PLOT_PROPS.max_grid_height,
            gridIndex: selectedSamples.length,
            show: false,
        });
    }
    return axes;
}

/**
 * Get max y-axis value
 * @param {WgsCovPlotDB} db
 * @returns {number}
 */
function getYAxisMax(db) {
    let yAxisMax = 0;
    if (db.segment_virus === true) {
        for (let segment of db.selectedSegments) {
            for (let sample of db.selectedSamples) {
                let maxDepth = db.mosdepth_info[sample][segment].max_depth;
                if (yAxisMax <= maxDepth) {
                    yAxisMax = maxDepth;
                }
            }
        }
    } else {
        for (let sample of db.selectedSamples) {
            let maxDepth = db.mosdepth_info[sample].max_depth;
            if (yAxisMax <= maxDepth) {
                yAxisMax = maxDepth;
            }
        }
    }
    return yAxisMax;
}

export {getXAxes, getYAxes, getYAxisMax};