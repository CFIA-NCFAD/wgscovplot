import {FLU_SEGMENT_COLOURS} from "../../util";
import {ECFeature, MaxSegmentLength, SegmentCoords, WgsCovPlotDB} from "../../db";
import {find, get, isNil, omitBy} from "lodash";


export function getCustomXAxisLabel(value: number, segments: string[], segCoords: SegmentCoords): string {
  for (const segment of segments) {
    const {start, end} = segCoords[segment];
    if (value >= start && value <= end) {
      const pos = value - start + 1
      return `${segment}:${pos.toLocaleString()}`;
    }
  }
  return ""
}

export const getMaxSegmentLength = (db: WgsCovPlotDB) => {
  const out: MaxSegmentLength = {};
  const selectedSegments = db.chartOptions.selectedSegments;
  if (isNil(selectedSegments)) {
    return out;
  }
  for (const segment of selectedSegments) {
    let maxLength = 0;
    for (const sample of db.chartOptions.selectedSamples) {
      const segmentsRefSeq = db.segments_ref_seq;
      const n = get(segmentsRefSeq, [sample, segment, "length"], 0);
      if (n > maxLength) {
        maxLength = n;
      }
    }
    out[segment] = maxLength;
  }
  /*Sample may miss some segments in analysis, remove those segment
  Handel the case when user select only samples that miss some segments in analysis
   */
  return omitBy(out, segLength => segLength === 0);
}

export const getSegmentCoords = (db: WgsCovPlotDB): SegmentCoords => {
  const out: SegmentCoords = {};
  let prev = {
    maxLength: 0,
    start: 0,
    end: 0,
  }
  const maxSegmentLength = db.maxSegmentLength;
  if (isNil(maxSegmentLength) || isNil(db.chartOptions.selectedSegments)) {
    return out;
  }
  for (const segment of db.chartOptions.selectedSegments) {
    const maxLength = maxSegmentLength[segment];
    const tempCoords = {
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
  const segCoords = db.segCoords;
  const selectedSegments = db.chartOptions.selectedSegments;
  if (isNil(segCoords) || isNil(selectedSegments) || selectedSegments.length === 0) {
    return '';
  }
  for (const segment of selectedSegments) {
    const {start, end} = segCoords[segment];
    if (position >= start && position <= end) {
      return segment
    }
  }
  return '';
}

export const getFluGeneFeature = (db: WgsCovPlotDB): ECFeature[] => {
  const geneFeature: ECFeature[] = [];
  let i = 0
  const segmentLength = db.maxSegmentLength;
  const segCoords = db.segCoords;
  if (isNil(segmentLength) || isNil(segCoords) || isNil(db.chartOptions.selectedSegments)) {
    return geneFeature;
  }
  for (const segment of db.chartOptions.selectedSegments) {
    const {start, end} = segCoords[segment];
    let color = get(FLU_SEGMENT_COLOURS, segment, "#000");
    if (!isNil(db.echart_features)) {
      const feature = find(db.echart_features, ["name", segment]);
      if (!isNil(feature)) {
        color = get(feature, ["itemStyle", "color"], get(FLU_SEGMENT_COLOURS, segment, "#000"));
      }
    }
    geneFeature.push({
      name: `${segment}`,
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
        "color": color
      }
    });
    i += 1;
  }
  return geneFeature;
}