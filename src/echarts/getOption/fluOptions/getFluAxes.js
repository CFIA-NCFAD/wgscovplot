
function getCustomXAxisLabel (value, segments, segmentsRange){
    let segment;
    let pos;
    for (let i=0; i < segmentsRange.length; i ++){
        if (value >= segmentsRange[i][0] && value <= segmentsRange[i][1]){
            pos = value - segmentsRange[i][0] + 1;
            segment = segments[i];
        }
    }
    return segment + ' - ' + pos.toLocaleString();
}

/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples name
 * @param {number} xAxisMax - Max value is set for x Axis
 * @param {boolean} isShowXAxisLabel - whether to show X Axis
 * @returns {Array<Object>}
 */
function getFluXAxes(samples, segments, xAxisMax, isShowXAxisLabel, segmentsRange) {
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
                formatter: function (value){
                    return getCustomXAxisLabel(value, segments, segmentsRange);
                }
            }
        });
    }
    // For gene feature plot
    axes.push({
        type: "value",
        gridIndex: samples.length,
        min: 1,
        max: xAxisMax,
        axisLabel: {
            interval: "auto",
            formatter: function (value){
                return getCustomXAxisLabel(value, segments, segmentsRange);
            }
        },
    });
    return axes;
}

export {getFluXAxes}