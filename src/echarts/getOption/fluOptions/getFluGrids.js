/**
 *  Define grid for the whole charts
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of segments
 * @returns {Array<Object>}
 */
function getFluGrids(samples, segments) {
    let grids = [];
    let padTop = 4.0; // Percentage
    let horizontalOffset = 4.0;
    let verticalOffset = 6.0;
    let padLeft = 4.0;
    let plotHeight = 90.0;
    let samplePlotHeight = (plotHeight - padTop)/samples.length - verticalOffset;
    let samplePlotWidth = (100.0 - padLeft)/segments.length - horizontalOffset;
    for (let i = 0; i < samples.length; i++){
        for (let j = 0; j < segments.length; j++){
            grids.push({
                show: true,
                top: padTop + i * (samplePlotHeight + verticalOffset) + '%',
                height: samplePlotHeight + '%',
                left: padLeft + j * (samplePlotWidth + horizontalOffset) + '%',
                width: samplePlotWidth + '%',
                borderColor: 'red'
            });
        }
    }
    console.log(grids)
    console.log(chart.getWidth())
    return grids;
}

export {getFluGrids};