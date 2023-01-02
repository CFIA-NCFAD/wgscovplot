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

export {
    getDataZoom
};