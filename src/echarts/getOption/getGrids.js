import {geneFeaturePlotProperties} from "../../util";

/**
 *  Define grid for the whole charts
 * @param {Array<string>} samples - An array of samples name
 * @param {boolean} geneFeature - whether to plot gene feature or not ("True" or "False")
 * @param {boolean} amplicon - whether to plot amplicon feature or not ("True" or "False")
 * @returns {Array<Object>}
 */
function getGrids(samples, geneFeature, amplicon) {
    var padTop = 4;
    var heightOffset = 6;
    var grids = [];
    var n = (amplicon || geneFeature) ? samples.length + 1 : samples.length;
    // Each subplot has same height, the height can be adjusted in control menu
    for (var idx = 0; idx < n; idx++){
        grids.push({
            show: true,
            height: (1/(n+1)) * 100 - heightOffset + "%",
            top: (idx/(n+1)) * 100 + padTop + "%",
            left: "8%",
            right: "8%",
        });
    }
    return grids;
}

export {getGrids};