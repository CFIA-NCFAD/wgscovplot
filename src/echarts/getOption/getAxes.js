/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples name
 * @param {number} xAxisMax - Max value is set for x Axis
 * @returns {Array<Dict[]>}
 */
function getXAxes(samples, xAxisMax) {
    var axes = [];
    for (var [i, sample] of samples.entries()) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                interval: "auto",
            },
        });
    }
    if (amplicon === "True" || geneFeature === "True"){
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
 * @returns {Array<Dict[]>}
 */
function getYAxes(samples, scaleType, yMax) {
    var axes = [];
    for (var [i, sample] of samples.entries()) {
        axes.push({
            type: scaleType,
            gridIndex: i,
            name: sample,
            nameTextStyle: {
                fontStyle: "normal",
                fontWeight: "bolder",
            },
            nameLocation: "end",
            min: scaleType === 'log' ? 1 : 0,
            max: yMax,
            minorSplitLine: {
                show: true,
            },
        });
    }
    if (amplicon === "True" || geneFeature === "True"){
        axes.push({
            max: geneFeatureProperties["max_grid_height"],
            gridIndex: samples.length,
            show: false,
        });
    }
    return axes;
}

export {getXAxes, getYAxes};