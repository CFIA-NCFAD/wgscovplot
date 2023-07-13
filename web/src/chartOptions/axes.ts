import {isNil} from "lodash";
import {SampleSegmentMosdepthInfo, SegmentCoords, WgsCovPlotDB} from "../db";


interface SegmentMaxLength {
  [key: string]: number;
}

interface SegmentInfo {
  [key: string]: {
    start: number;
    end: number;
    maxLength: number;
  }
}

const getMaxSegmentLength = (db: WgsCovPlotDB): SegmentMaxLength => {
  let selectedSegments = db.chartOptions.selectedSegments;
  if (isNil(selectedSegments)) {
    selectedSegments = db.segments;
  }
  if (selectedSegments === undefined) {
    return {};
  }
  let selectedSamples = db.chartOptions.selectedSamples;
  let out: SegmentMaxLength = {}
  const depthInfo = db.mosdepth_info as SampleSegmentMosdepthInfo;
  for (let segment of selectedSegments) {
    let maxLength = 0;

    for (let sample of selectedSamples) {
      if (depthInfo[sample] === undefined) {
        continue;
      }
      if (depthInfo[sample][segment] === undefined) {
        continue;
      }
      let length = depthInfo[sample][segment].ref_seq_length;
      maxLength = Math.max(maxLength, length)
    }
    out[segment] = maxLength;
  }
  return out
}


export const getSegmentCoords = (db: WgsCovPlotDB): SegmentInfo => {
  let selectedSegments = db.chartOptions.selectedSegments;
  if (isNil(selectedSegments)) {
    selectedSegments = db.segments;
  }
  if (selectedSegments === undefined) {
    return {};
  }
  let out: SegmentInfo = {};
  let prev = {
    maxLength: 0,
    start: 0,
    end: 0,
  }
  let maxSegmentLength = getMaxSegmentLength(db);
  for (let segment of selectedSegments) {
    let maxLength = maxSegmentLength[segment];
    let segCoords = {
      maxLength,
      start: prev.end + 1,
      end: prev.end + maxLength,
    };
    out[segment] = segCoords;
    prev = segCoords;
  }
  return out;
}

/**
 * Custom xAxis label (for segment virus)
 * @param {number} value - xAxis value
 * @param {Array<string>} segments - An array of segments names
 * @param {SegmentCoords} segCoords
 * @returns {string}
 */
export function getCustomXAxisLabel(value: number, segments: string[], segCoords: SegmentCoords): string {
  for (let segment of segments) {
    let {start, end} = segCoords[segment];
    if (value >= start && value <= end) {
      let pos = value - start + 1
      return `${segment}:${pos.toLocaleString()}`;
    }
  }
  return ""
}
