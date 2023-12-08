import {WgsCovPlotDB} from "../db";

export const getGrids = (db: WgsCovPlotDB) => {
  console.log("Trigger getGrids")
  let featureHeight: number = 0.0;
  if (db.chartOptions.showFeatures && db.show_amplicons && db.show_genes) {
    featureHeight = 15.0;
  } else if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    featureHeight = 6.0;
  } else {
    // no features subplot shown
    featureHeight = -5.0;
  }
  // TODO: find out what doubleStrand is
  featureHeight = (db.doubleStrand && featureHeight > 0) ? (featureHeight + 6.0) : featureHeight;
  featureHeight *= (db.chartOptions.featurePlotHeightScaling / 100);
  let padTop = db.chartOptions.padTop;
  let plotHeight = (db.chartOptions.showDataZoomSlider) ? 90 : 95;
  let subPlotHeight = plotHeight - featureHeight;
  let grids = [];
  let nSamples = db.chartOptions.selectedSamples.length;
  let heightOffset = db.chartOptions.heightOffset;
  let sampleHeight = (subPlotHeight - padTop) / nSamples - heightOffset;
  for (let idx = 0; idx < nSamples; idx++) {
    grids.push({
      show: true,
      height: sampleHeight + "%",
      top: ((sampleHeight + heightOffset) * idx + padTop) + "%",
      left: db.chartOptions.leftMargin + "%",
      right: db.chartOptions.rightMargin + "%",
    });
  }
  if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    grids.push({
      show: false,
      height: featureHeight + "%",
      top: ((sampleHeight + heightOffset) * nSamples + padTop) + "%",
      left: db.chartOptions.leftMargin + "%",
      right: db.chartOptions.rightMargin + "%",
    });
  }
  return grids;
}
