import {geneFeaturePlotProperties} from "../../util";

/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples name
 * @param {number} xAxisMax - Max value is set for x Axis
 * @param {boolean} geneFeature - whether to plot gene feature or not (true for false)
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @param {boolean} isShowXAxisLabel - whether to show X Axis
 * @returns {Array<Object>}
 */
function getXAxes(samples, xAxisMax, geneFeature, amplicon, isShowXAxisLabel) {
    let axes = [];
    for (let i = 0; i < samples.length; i++) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                show:isShowXAxisLabel,
                interval: "auto",
            }
        });
    }
    if (amplicon || geneFeature){
        axes.push({
            type: "value",
            gridIndex: samples.length,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                interval: "auto",
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
            min: scaleType === 'log' ? 1 : 0,
            max: yMax,
            minorSplitLine: {
                show: true,
            },
        });
    }
    if (amplicon || geneFeature){
        axes.push({
            max: geneFeaturePlotProperties.max_grid_height,
            gridIndex: samples.length,
            show: false,
        });
    }
    return axes;
}

export {getXAxes, getYAxes};