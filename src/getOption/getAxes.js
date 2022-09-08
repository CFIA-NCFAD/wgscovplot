import {geneFeaturePlotProperties} from "../util";

/**
 * Custom xAxis label (for segment virus)
 * @param {number} xAxisValue - xAxis value
 * @param {Array<string>} segments - An array of segments names
 * @param {Array<Array<number>>} segmentsInterval - An array of segment start, end
 * @returns {string}
 */
function getCustomXAxisLabel(xAxisValue, segments, segmentsInterval) {
    let segment = "";
    let pos = 1;
    for (let i = 0; i < segmentsInterval.length; i++) {
        if (xAxisValue >= segmentsInterval[i][0] && xAxisValue <= segmentsInterval[i][1]) {
            pos = xAxisValue - segmentsInterval[i][0] + 1;
            segment = segments[i];
        }
    }
    return segment + " - " + pos.toLocaleString();
}

/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names (param for segment virus)
 * @param {Array<Array<number>>} segmentsInterval - An array of segment start, end (param for segment virus)
 * @param {number} xAxisMax - Max value is set for x Axis
 * @param {boolean} geneFeature - whether to plot gene feature or not (true for false)
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @param {boolean} showXAxisLabel - whether to show X Axis
 * @returns {Array<Object>}
 */
function getXAxes(samples, segments, segmentsInterval, xAxisMax,
                  geneFeature, amplicon, showXAxisLabel) {
    let axes = [];
    for (let i = 0; i < samples.length; i++) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                show: showXAxisLabel,
                interval: "auto",
                formatter: (segments.length > 0) ? function (value) {
                    return getCustomXAxisLabel(value, segments, segmentsInterval);
                } : {}
            }
        });
    }
    if (amplicon || geneFeature) {
        axes.push({
            type: "value",
            gridIndex: samples.length,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                interval: "auto",
                formatter: (segments.length > 0) ? function (value) {
                    return getCustomXAxisLabel(value, segments, segmentsInterval);
                } : {}
            },
        });
    }
    return axes;
}

/**
 * Define options for Y axis
 * @param {Array<string>} samples - An array of samples name
 * @param {string} scaleType - scale for Y Axis, either value or log
 * @param {number} yMax - max value is set for y Axis
 * @param {boolean} geneFeature - whether to plot gene feature or not (true for false)
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @returns {Array<Object>}
 */
function getYAxes(samples, scaleType, yMax, geneFeature, amplicon) {
    let axes = [];
    for (let [i, sample] of samples.entries()) {
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
    if (amplicon || geneFeature) {
        axes.push({
            max: geneFeaturePlotProperties.max_grid_height,
            gridIndex: samples.length,
            show: false,
        });
    }
    return axes;
}

export {getXAxes, getYAxes};