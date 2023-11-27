import {FLU_SEGMENT_COLOURS} from "../../util";
import {max, sum, values} from "lodash";
import {ECFeature, MaxSegmentLength, SampleSegmentDepths, SegmentCoords, WgsCovPlotDB} from "../../db";


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

export const getMaxSegmentLength = (db: WgsCovPlotDB) => {
  let out: MaxSegmentLength = {};
  for (let segment of db.chartOptions.selectedSegments) {
    let maxLength = 0;
    for (let sample of db.chartOptions.selectedSamples) {
      // @ts-ignore
      let length = db.mosdepth_info[sample][segment].ref_seq_length;
      if (maxLength <= length) {
        maxLength = length;
      }
    }
    out[segment] = maxLength;
  }
  return out
}


export const getSegmentCoords = (db: WgsCovPlotDB): SegmentCoords => {
  let out: SegmentCoords = {};
  let prev = {
    maxLength: 0,
    start: 0,
    end: 0,
  }
  for (let segment of db.chartOptions.selectedSegments) {
    let maxLength = db.maxSegmentLength[segment];
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

export const whichSegment = (position: number, db: WgsCovPlotDB) : string => {
    for (let segment of Object.keys(db.segCoords)) {
        let coords = db.segCoords[segment];
        if (position >= coords.start && position <= coords.end) {
            return segment
        }
    }
    return '';
}


export const getFluGeneFeature = (db: WgsCovPlotDB) : ECFeature[] => {
  let geneFeature: ECFeature[] = [];
  let i = 0
  for (let segment of db.chartOptions.selectedSegments) {
    let {start, end} = db.segCoords[segment];
    geneFeature.push({
      name: `Segment ${segment}`,
      value: {
        idx: i,
        start,
        end,
        level: 0,
        strand: 1,
        rotate: 0.0,
        type: "segment",
      },
      itemStyle: {
        // @ts-ignore
        "color": FLU_SEGMENT_COLOURS[segment]
      }
    });
    i += 1;
  }
  return geneFeature;
}