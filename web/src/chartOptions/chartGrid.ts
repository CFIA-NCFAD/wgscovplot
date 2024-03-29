import {WgsCovPlotDB} from "../db";

export const getGrids = (db: WgsCovPlotDB) => {
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
  const padTop = db.chartOptions.padTop;
  const plotHeight = (db.chartOptions.showDataZoomSlider) ? 90 : 95;
  const subPlotHeight = plotHeight - featureHeight;
  const grids = [];
  const nSamples = db.chartOptions.selectedSamples.length;
  const heightOffset = db.chartOptions.heightOffset;
  const sampleHeight = (subPlotHeight - padTop) / nSamples - heightOffset;
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
