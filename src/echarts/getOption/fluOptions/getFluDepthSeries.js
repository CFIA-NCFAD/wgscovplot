
function getFluDepthSeries(samples, segments) {
    let depthSeries = [];
    for (let i = 0; i < samples.length * segments.length; i++) {
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

export {getFluDepthSeries}