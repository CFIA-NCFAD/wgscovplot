/**
 * Get data zoom option
 * @param {Array<string>} samples - An array of samples names
 * @returns {Array<Object>}
 */
function getDataZoom(samples) {
    let dataZoom;
    dataZoom = [
        {
            type: "inside",
            filterMode: "none",
            xAxisIndex: [...Array(samples.length + 1).keys()],
            zoomLock: false,
        },
        {
            show: true,
            filterMode: "none",
            xAxisIndex: [...Array(samples.length + 1).keys()],
            type: "slider",
            zoomLock: false,
            showDataShadow: false
        },
    ];
    return dataZoom;
}

export {getDataZoom};