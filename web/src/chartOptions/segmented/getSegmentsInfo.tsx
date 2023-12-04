import {FLU_SEGMENT_COLOURS} from "../../util";
import {ECFeature, MaxSegmentLength, SegmentCoords, WgsCovPlotDB} from "../../db";
import {omitBy} from "lodash";


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
      let length = db.segments_ref_seq[sample][segment].length;
      if (maxLength <= length) {
        maxLength = length;
      }
    }
    out[segment] = maxLength;
  }
  /*Sample may miss some segments in analysis, remove those segment
  Handel the case when user select only samples that miss some segments in analysis
   */
  let filteredResults: MaxSegmentLength = omitBy(out, segLength => segLength === 0)
  return filteredResults
}

export const getSegmentCoords = (db: WgsCovPlotDB): SegmentCoords => {
  let out: SegmentCoords = {};
  let prev = {
    maxLength: 0,
    start: 0,
    end: 0,
  }
  for (let segment of Object.keys(db.maxSegmentLength)) {
    let maxLength = db.maxSegmentLength[segment];
    let tempCoords = {
      maxLength,
      start: prev.end + 1,
      end: prev.end + maxLength,
    };
    out[segment] = tempCoords;
    prev = tempCoords;
  }
  return out;
}

export const whichSegment = (position: number, db: WgsCovPlotDB): string => {
  for (let segment of Object.keys(db.segCoords)) {
    let {start, end} = db.segCoords[segment];
    if (position >= start && position <= end) {
      return segment
    }
  }
  return '';
}

export const getFluGeneFeature = (db: WgsCovPlotDB): ECFeature[] => {
  let geneFeature: ECFeature[] = [];
  let i = 0
  for (let segment of Object.keys(db.maxSegmentLength)) {
    let {maxLength, start, end} = db.segCoords[segment];
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