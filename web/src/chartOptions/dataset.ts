import {ECFeature, MaxSegmentLength, SegmentCoords, WgsCovPlotDB} from "../db";
import {constant, isNil, sum, times, values} from "lodash";
import {getFluGeneFeature, getMaxSegmentLength, getSegmentCoords} from "./segmented/getSegmentsInfo";
import {setState, state} from "../state";
import {leaveSelect} from "echarts/types/src/util/states";

export const getDatasets = (db: WgsCovPlotDB) => {
  let datasets = [];
  if (!isNil(db.segments)) { // segmented virus
    /* For nf-flu, setup must be done before setting chart option:
    Set up gene feature, position array and segment coords
    */
    let maxSegmentLength: MaxSegmentLength = getMaxSegmentLength(state);
    setState({"maxSegmentLength": maxSegmentLength}) // disable shallow merging and fully replace
    let positions: number[] = [...Array(sum(values(maxSegmentLength)) + 1).keys()];
    positions.shift()
    setState("positions", positions);
    let segCoords: SegmentCoords = getSegmentCoords(state);
    setState({"segCoords": segCoords}); // disable shallow merging and fully replace
    let echart_features: ECFeature[] = getFluGeneFeature(state);
    setState("echart_features", echart_features)
    //////////////// Done Setup //////////////////////////////

    for (let sample of db.chartOptions.selectedSamples) {
      let depthArray: number [] = [];
      for (let segment of Object.keys(segCoords)) {
        // @ts-ignore
        let ds: number[] = isNil(db.depths[sample][segment]) ? [] : db.depths[sample][segment];
        let coords = db.segCoords[segment];
        if (ds.length < coords.maxLength) {
          let padding = times(coords.maxLength - ds.length, constant(1E-7));
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
  } else { // non-segmented virus
    for (let sample of db.chartOptions.selectedSamples) {
      datasets.push({
        dimensions: [
          {name: "depth", type: "float"},
          {name: "position", type: "int"},
        ],
        source: {
          position: db.positions,
          depth: db.depths[sample],
        },
      });
    }
  }
  return datasets;
}