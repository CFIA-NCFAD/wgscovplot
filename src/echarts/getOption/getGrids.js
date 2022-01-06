import {geneFeaturePlotProperties} from "../../util";

/**
 *  Define grid for the whole charts
 * @param {Array<string>} samples - An array of samples name
 * @param {boolean} geneFeature - whether to plot gene feature or not ("True" or "False")
 * @param {boolean} amplicon - whether to plot amplicon feature or not ("True" or "False")
 * @returns {Array<Object>}
 */
function getGrids(samples, geneFeature, amplicon) {
    var n = samples.length + 1;
    var lastHeight;
    var lastTop;
    var grids = Object.keys(samples).map(function (sample) {
        lastHeight = (1 / n) * 100 - 6;
        if (n === 2) {
            // Only 1 sample (1 sample + gene feature plot)
            lastHeight = 70;
            return {
                show: true,
                height: "70%", // plot display in nearly full scale
            };
        }
        return {
            show: true,
            height: (1 / n) * 100 - 6 + "%",
        };
    });
    grids.forEach(function (grid, idx) {
        var padTop = 4;
        lastTop = (idx / n) * 100 + padTop;
        grid.top = (idx / n) * 100 + padTop + "%";
        grid.left = "8%";
        grid.right = "8%";
    });
    if (amplicon || geneFeature){
        grids.push({
            show: true,
            height: geneFeaturePlotProperties["grid_height"],
            top: lastHeight + lastTop + 3 + "%",
            left: "8%",
            right: "8%",
        });
    }
    return grids;
}

export {getGrids};