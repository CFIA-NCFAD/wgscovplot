import {WgsCovPlotDB} from "../db";
import {constant, get, isNil, times} from "lodash";

export const getDatasets = (db: WgsCovPlotDB) => {
  const datasets = [];
  const depths = db.depths;
  if (!isNil(db.segments) && !isNil(db.chartOptions.selectedSegments)) {
    // segmented virus
    for (const sample of db.chartOptions.selectedSamples) {
      let depthArray: number [] = [];
      for (const segment of db.chartOptions.selectedSegments) {
        let ds: number[] = get(depths, [sample, segment], []);
        const maxLength = get(db.maxSegmentLength, segment, 0);
        if (ds.length < maxLength) {
          const padding = times(maxLength - ds.length, constant(1E-7));
          ds = [...ds, ...padding]
        }
        depthArray = [...depthArray, ...ds];
      }
      datasets.push({
        dimensions: [
          {name: "depth", type: "float"},
          {name: "position", type: "int"},
        ],
        source: {
          position: db.positions,
          depth: depthArray,
        },
      });
    }
  } else {
    // non-segmented virus
    for (const sample of db.chartOptions.selectedSamples) {
      datasets.push({
        dimensions: [
          {name: "depth", type: "float"},
          {name: "position", type: "int"},
        ],
        source: {
          position: db.positions,
          depth: depths[sample],
        },
      });
    }
  }
  return datasets;
}