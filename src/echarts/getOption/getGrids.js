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
    var verticalRatio = 0.5;
    var grids = [];
    var n = (amplicon || geneFeature) ? samples.length + 1 : samples.length;
    // Each subplot has same height, the height can be adjusted in control menu
    for (var idx = 0; idx < n; idx++){
        grids.push({
            show: true,
            height: (1/(n+verticalRatio)) * 100 - heightOffset + "%", // vertical space
            top: (idx/(n+verticalRatio)) * 100 + padTop + "%", // vertical space
            left: "4%", // horizontal space
            right: "4%", // horizontal space
        });
    }
    return grids;
}

export {getGrids};