/**
 * Custom xAxis label
 * @param {number} xAxisValue - xAxis value
 * @param {Array<string>} segments - An array of segments names
 * @param {Array<Array<number>>} segmentsRange - An array of segment start, end
 * @returns {string}
 */
function getCustomXAxisLabel(xAxisValue, segments, segmentsRange) {
    let segment = '';
    let pos = 1;
    for (let i = 0; i < segmentsRange.length; i++) {
        if (xAxisValue >= segmentsRange[i][0] && xAxisValue <= segmentsRange[i][1]) {
            pos = xAxisValue - segmentsRange[i][0] + 1;
            segment = segments[i];
        }
    }
    return segment + ' - ' + pos.toLocaleString();
}

/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {number} xAxisMax - Max value is set for xAxis
 * @param {boolean} showXAxisLabel - whether to show X Axis label
 * @param {Array<Array<number>>} segmentsRange - An array of segment start, end
 * @returns {Array<Object>}
 */
function getFluXAxes(samples, segments, xAxisMax,
                     showXAxisLabel, segmentsRange) {
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
                formatter: function (value) {
                    return getCustomXAxisLabel(value, segments, segmentsRange);
                }
            }
        });
    }
    axes.push({ // For gene feature plot
        type: "value",
        gridIndex: samples.length,
        min: 1,
        max: xAxisMax,
        axisLabel: {
            interval: "auto",
            formatter: function (value) {
                return getCustomXAxisLabel(value, segments, segmentsRange);
            }
        },
    });
    return axes;
}

export {getFluXAxes};