/**
 * Get dataset information
 * @param {Array<Array<number>>} depths - Array of depths
 * @returns {Array<Object>}
 */
function getFluDatasets(depths) {
    let datasets = [];
    for (let [i, depthArray] of depths.entries()) {
        let position = [...Array(depthArray.length + 1).keys()]
        position.shift();
        datasets.push({
            dimensions: [
                {name: "depth", type: "float"},
                {name: "position", type: "int"},
            ],
            source: {
                position: position,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

export {getFluDatasets};