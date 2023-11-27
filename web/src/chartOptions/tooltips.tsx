import {get, isNil, find, uniqBy, flatMap, keys} from "lodash";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {SampleDepths, VariantCall, WgsCovPlotDB} from "../db";
import {unwrap} from "solid-js/store";


export const getVariantComparison = (db: WgsCovPlotDB, position: number): Array<Array<string>> => {
  if (isNil(db.variants)) {
    return [];
  }
  let rows: any[][] = [];
  let variantArr: VariantCall[] = [];
  for (let sample of db.chartOptions.selectedSamples) {
    let variant = find(db.variants[sample], {POS: position}, 0) as VariantCall;
    if (!isNil(variant)) {
      variantArr.push(variant);
    } else {
      variantArr.push({sample, "POS": position.toString()});
    }
  }
  // @ts-ignore
  let unionKeys = uniqBy(flatMap(variantArr, keys));
  let depths: SampleDepths = unwrap(db.depths) as SampleDepths;
  unionKeys.push("Coverage Depth");
  unionKeys.forEach((key: string) => {
    let row = [];
    row.push(key);
    if (key === "Coverage Depth") {
      for (let sample of db.chartOptions.selectedSamples) {
        row.push(depths[sample][position - 1].toLocaleString());
      }
    } else {
      for (let variant of variantArr) {
        let val = get(variant, key, "");
        row.push(val);
      }
    }
    rows.push(...[row]);
  });
  return rows;
}

export const getCoverageStatComparison = (db: WgsCovPlotDB, start: number, end: number, position: number) => {
  const selectedSamples = db.chartOptions.selectedSamples;
  const depths = unwrap(db.depths) as SampleDepths;
  const low_coverage_threshold = db.chartOptions.low_coverage_threshold;
  let rows = [];
  let tableHeader = [
    "Sample",
    `Depth at position ${position.toLocaleString()}`,
    "Range",
    "Mean Coverage (X)",
    "Median Coverage (X)",
    `Genome Coverage (>=${low_coverage_threshold}X) (%)`,
  ];
  rows.push(...[tableHeader]);
  for (let sample of selectedSamples) {
    let sampleDepths = depths[sample];
    let meanCov = meanCoverage(sampleDepths, start, end).toFixed(2);
    let medianCov = medianCoverage(sampleDepths, start, end).toFixed(2);
    let genomeCov = genomeCoverage(sampleDepths, start, end, low_coverage_threshold).toFixed(2);
    let coverageDepth = sampleDepths[position - 1];
    let row = [
      sample,
      coverageDepth.toLocaleString(),
      `${start.toLocaleString()}-${end.toLocaleString()}`,
      meanCov,
      medianCov,
      genomeCov,
    ];
    rows.push(...[row]);
  }
  return rows;
}

export const getSegmentVariantComparison = (db: WgsCovPlotDB, sample: string, segment: string, position: number): Array<Array<string>> => {
  let variantArr: any = [];
  db.chartOptions.selectedSamples.forEach(sample => {
    // @ts-ignore
    let variantCall = find(db.variants[sample][segment],
      {"POS": position.toString()}, 0);
    if (!isNil(variantCall)) {
      variantArr.push(variantCall);
    } else {
      variantArr.push({
        "Sample": sample,
        "POS": position,
        "REF_ID": db.segments_ref_id[sample][segment],
        "Segment": segment,
        // @ts-ignore
        "Segment Length": db.depths[sample][segment].length,
        "REF_SEQ": db.segments_ref_seq[sample][segment][position - 1],
      });
    }
  });
  // @ts-ignore
  let unionKeys = uniqBy(flatMap(variantArr, keys));
  unionKeys.push("Coverage Depth"); // Add Coverage Depth row
  let rows: string[][] = [];
  unionKeys.forEach(key => {
    let row = [];
    row.push(key);
    if (key === "Coverage Depth") {
      db.chartOptions.selectedSamples.forEach(sample => {
        // @ts-ignore
        const segDepths = db.depths[sample][segment];
        let rowValue;
        if (segDepths.length === 0) {
          rowValue = `No result reported for segment ${segment}`;
        } else {
          if (position <= segDepths.length) {
            rowValue = segDepths[position - 1].toLocaleString();
          } else {
            //Out of range when segment length < padding array
            rowValue = `No sequence at this position. Reference sequence 
                            ${db.segments_ref_id[sample][segment]} is only 
                            ${segDepths.length} bp`;
          }
        }
        row.push(rowValue);
      });
    } else {
      variantArr.forEach((element: { [x: string]: any; }) => {
        let variant = element[key];
        if (!isNil(variant)) {
          row.push(variant);
        } else {
          row.push(""); // No information
        }
      });
    }
    rows.push(...[row]);
  });
  return rows;
}

export const getSegmentCoverageStatComparison = (db: WgsCovPlotDB, segment: string, position: number) => {
  let rows = [];
  let tableHeader = [
    "Sample",
    `Depth at position ${position}`,
    "Segment",
    "Range",
    "Segment Length",
    "Mean Coverage (X)",
    "Median Coverage (X)",
    `Genome Coverage (>=${db.chartOptions.low_coverage_threshold}X) (%)`
  ];
  rows.push(...[tableHeader]);
  db.chartOptions.selectedSamples.forEach((sample) => {
    // @ts-ignore
    let depthArr = db.depths[sample][segment];
    let meanCov = meanCoverage(depthArr, 1, depthArr.length).toFixed(2);
    let medianCov = medianCoverage(depthArr, 1, depthArr.length).toFixed(2);
    let genomeCov = genomeCoverage(depthArr, 1, depthArr.length, db.chartOptions.low_coverage_threshold).toFixed(2);
    let coverageDepth;
    if (position <= depthArr.length) {
      coverageDepth = depthArr[position - 1].toLocaleString();
    } else if (depthArr.length === 0) {
      coverageDepth = `No result reported for segment ${segment}`;
    } else {
      coverageDepth = "No sequence at this position";
    }
    let row = [
      sample,
      coverageDepth,
      segment,
      `1-${depthArr.length}`,
      depthArr.length,
      meanCov,
      medianCov,
      genomeCov
    ];
    rows.push(row);
  });
  return rows;
}