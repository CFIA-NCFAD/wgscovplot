/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples name
 * @returns {Array<Object>}
 */
function getDepthSeries(samples) {
    var depthSeries = [];
    for (var [i, sample] of samples.entries()) {
        depthSeries.push({
            type: "line",
            xAxisIndex: i,
            yAxisIndex: i,
            areaStyle: {
                color: "#666",
            },
            encode: {
                x: "position",
                y: "depth",
            },
            symbol: "none",
            datasetIndex: i,
            lineStyle: {
                color: "#666",
                opacity: 0,
            },
            large: true,
        });
    }
    return depthSeries;
}

export {getDepthSeries}