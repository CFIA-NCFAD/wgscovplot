/**
 * Get ECharts Grid options objects given selected samples (and segments) and if gene/amplicon/segment features are being shown
 * @param {WgsCovPlotDB} db
 * @returns {Object[]}
 */
function getGrids(db) {
    const {
        selectedSamples,
        show_genes,
        show_amplicons,
        doubleStrand,
        padTop = 4.0,
        heightOffset = 6.0,
    } = db;
    let featureHeight;
    // TODO: make magic numbers function optional params
    if (show_amplicons && show_genes) {
        featureHeight = 15.0;
    } else if (show_amplicons || show_genes) {
        featureHeight = 6.0;
    } else {
        featureHeight = -5.0; // make subplot full
    }
    // TODO: find out what doubleStrand is
    featureHeight = (doubleStrand && featureHeight > 0) ? (featureHeight + 6.0) : featureHeight;
    let subPlotHeight = 90.0 - featureHeight;
    let grids = [];
    let sampleHeight = (subPlotHeight - padTop) / selectedSamples.length - heightOffset;
    for (let idx = 0; idx < selectedSamples.length; idx++) {
        grids.push({
            show: true,
            height: sampleHeight.toFixed(1) + "%",
            top: ((sampleHeight + heightOffset) * idx + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    if (show_amplicons || show_genes) {
        grids.push({
            show: false,
            height: featureHeight.toFixed(1) + "%",
            top: ((sampleHeight + heightOffset) * selectedSamples.length + padTop).toFixed(1) + "%",
            left: "4.0%",
            right: "4.0%",
        });
    }
    return grids;
}

export {
    getGrids
};