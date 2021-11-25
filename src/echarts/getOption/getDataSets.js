/**
 * Get dataset information
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<number>} positions - an array of genome positions which represent in X Axis
 * @returns {Array<Dict[]>}
 */
function getDatasets(depths, positionArray) {
    var datasets = [];
    for (var [i, depthArray] of depths.entries()) {
        datasets.push({
            dimensions: [
                {name: "depth", type: "float"},
                {name: "position", type: "int"},
            ],
            source: {
                position: positionArray,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

export {getDatasets};