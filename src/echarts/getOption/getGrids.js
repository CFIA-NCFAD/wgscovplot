import {geneFeaturePlotProperties} from "../../util";

/**
 *  Define grid for the whole charts
 * @param {Array<string>} samples - An array of samples name
 * @param {boolean} geneFeature - whether to plot gene feature or not
 * @param {boolean} amplicon - whether to plot amplicon feature or not
 * @param {boolean} doubleStrand - whether to plot amplicon feature or not
 * @returns {Array<Object>}
 */
function getGrids(samples, geneFeature, amplicon, doubleStrand) {
    var padTop = 4.0; // Percentage
    var heightOffset = 6.0;
    var featureHeight;
    var subPlotHeight;
    if (amplicon && geneFeature){
        featureHeight = 15.0;
    } else if (amplicon || geneFeature){
        featureHeight = 6.0;
    } else{
        featureHeight = -5.0; // make subplot full
    }
    featureHeight = (doubleStrand && featureHeight > 0) ? (featureHeight + 6.0) : featureHeight;
    subPlotHeight = 90.0 - featureHeight;
    var grids = [];
    var sampleHeight = (subPlotHeight - padTop) / samples.length - heightOffset;
    for (var idx = 0; idx < samples.length; idx++){
        grids.push({
            show: true,
            height: parseFloat(sampleHeight).toFixed(1) + "%",
            top: parseFloat((sampleHeight + heightOffset) * idx + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    if (amplicon || geneFeature) {
        grids.push({
            show: false,
            height: parseFloat(featureHeight).toFixed(1) + "%",
            top: parseFloat((sampleHeight + heightOffset) * samples.length + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    return grids;
}

export {getGrids};