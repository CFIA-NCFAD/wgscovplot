/**
 * Get dataset information
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<number>} positions - an array of genome positions which represent in X Axis
 * @returns {Array<Object>}
 */
function getDatasets(depths, positions) {
    let datasets = [];
    for (let [i, depthArray] of depths.entries()) {
        datasets.push({
            dimensions: [
                {name: "depth", type: "float"},
                {name: "position", type: "int"},
            ],
            source: {
                position: positions,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

export {getDatasets};