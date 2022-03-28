/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples name
 * @param {boolean} isNonVariantSites - whether to show tooltips for non variant sites
 * @returns {Array<Object>}
 */
function getDepthSeries(samples, isNonVariantSites) {
    let depthSeries = [];
    for (let i = 0; i < samples.length; i++) {
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
            tooltip:{
                trigger: isNonVariantSites ? "axis" : "none"
            },
            silent: true,
            large: true,
        });
    }
    return depthSeries;
}

export {getDepthSeries}