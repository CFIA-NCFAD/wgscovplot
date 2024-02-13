import {find, flatMap, get, isEmpty, isNil, keys, uniq} from "lodash";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {SampleDepths, SampleSegmentDepths, VariantCall, WgsCovPlotDB} from "../db";
import {unwrap} from "solid-js/store";
import {setState} from "../state";
import {whichSegment} from "./segmented/getSegmentsInfo";
import {getPrimerInfo} from "./segmented/pcr";


export function getVariantComparison(
  db: WgsCovPlotDB,
  position: number,
  samples: string[]
) {
  if (isEmpty(db.variants)) {
    return [];
  }
  const rows: (string | number)[][] = [];
  const variantArr: VariantCall[] = [];
  for (const sample of samples) {
    const variant = find(db.variants[sample], {POS: position}, 0) as VariantCall;
    const depth = get(unwrap(db.depths), [sample, position - 1], 0) as number;
    if (!isNil(variant)) {
      variantArr.push({
        ...variant,
        "Sample": sample,
        "Coverage Depth": depth.toLocaleString()
      });
    } else {
      variantArr.push({
        "Sample": sample,
        "POS": position,
        "ALT": "N/A",
        "REF": "N/A",
        "Coverage Depth": depth.toLocaleString(),
      });
    }
  }
  const allKeys = [
    "Sample",
    "POS",
    "REF",
    "ALT",
    "Coverage Depth",
    ...flatMap(variantArr, keys)
  ];
  const unionKeys = uniq(allKeys);
  for (const key of unionKeys) {
    if (key === "sample") {
      continue;
    }
    const row = [];
    row.push(key);
    for (const variant of variantArr) {
      row.push(get(variant, key, ""));
    }
    rows.push(row);
  }
  return rows;
}

export function getCoverageStatComparison(
  db: WgsCovPlotDB,
  start: number,
  end: number,
  position: number
) {
  const selectedSamples = db.chartOptions.selectedSamples;
  const depths = unwrap(db.depths) as SampleDepths;
  const low_coverage_threshold = db.chartOptions.low_coverage_threshold;
  const rows: string[][] = [];
  const headers = [
    "Sample",
    `Depth at position ${position.toLocaleString()}`,
    "Range",
    "Mean Coverage (X)",
    "Median Coverage (X)",
    `Genome Coverage (>=${low_coverage_threshold}X) (%)`,
  ];
  rows.push(headers);
  for (const sample of selectedSamples) {
    const sampleDepths = depths[sample];
    const meanCov = meanCoverage(sampleDepths, start, end).toFixed(2);
    const medianCov = medianCoverage(sampleDepths, start, end).toFixed(2);
    const genomeCov = genomeCoverage(sampleDepths, start, end, low_coverage_threshold).toFixed(2);
    const coverageDepth = sampleDepths[position - 1];
    const row = [
      sample,
      coverageDepth.toLocaleString(),
      `${start.toLocaleString()}-${end.toLocaleString()}`,
      meanCov,
      medianCov,
      genomeCov,
    ];
    rows.push(row);
  }
  return rows;
}

interface SegmentedVariantComparison {
  [key: string]: string | number;

  Sample: string;
  POS: number;
  Segment: string;
  "Segment Length": number;
  REF_ID: string;
  REF_SEQ: string;
}

function missingSegmentVariantInfo(sample: string, position: number, segment: string) {
  return {
    "Sample": sample,
    "POS": position.toString(),
    "Segment": segment,
    "Segment Length": "N/A",
    "REF_ID": "N/A",
    "REF_SEQ": "N/A",
    "Coverage Depth": "N/A",
    "ALT_SEQ": "N/A",
    "ALT_FREQ": "N/A",
  };
}

export function getSegmentVariantComparison(
  db: WgsCovPlotDB,
  samples: string[],
  segment: string,
  position: number,
) {
  const variantArr: (VariantCall | SegmentedVariantComparison)[] = [];
  for (const sample of samples) {
    const variants = get(db, ["variants", sample, segment]);
    if (isNil(variants)) {
      variantArr.push(missingSegmentVariantInfo(sample, position, segment));
      continue;
    }
    const variant = find(variants, {POS: position}, 0) as VariantCall;
    const ref_id = get(db, ["segments_ref_id", sample, segment], "");
    const ref_seq = get(db, ["segments_ref_seq", sample, segment, position - 1], "");
    const seq_depths = get(db, ["depths", sample, segment], []);
    const seg_len = seq_depths.length;
    if (!isNil(variant)) {
      variantArr.push({...variant, "Coverage Depth": seq_depths[position - 1]});
    } else {
      if (ref_id === "" || ref_seq === "" || seg_len === 0) {
        variantArr.push(missingSegmentVariantInfo(sample, position, segment));
        continue;
      }
      variantArr.push({
        "Sample": sample,
        "POS": position,
        "REF_ID": ref_id,
        "Segment": segment,
        "Segment Length": seg_len,
        "REF_SEQ": ref_seq,
        "ALT_SEQ": "N/A",
        "ALT_FREQ": "N/A",
        "Coverage Depth": seq_depths[position - 1],
      });
    }
  }
  const ordered = [
    "Sample",
    "Segment",
    "POS",
    "Segment Length",
    "Coverage Depth",
    "REF_ID",
    "REF_SEQ",
    "ALT_SEQ",
    "ALT_FREQ"
  ];
  const allKeys = [...ordered, ...flatMap(variantArr, keys)];
  const unionKeys = uniq(allKeys);
  const rows: (string | number)[][] = [];
  for (const key of unionKeys) {
    const row: (string | number)[] = [];
    row.push(key);
    for (const variant of variantArr) {
      row.push(get(variant, [key], ""));
    }
    rows.push(row);
  }
  return rows;
}

export function getSegmentCoverageStatComparison(
  db: WgsCovPlotDB,
  segment: string,
  position: number
) {
  const rows = [];
  const tableHeader = [
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
    const depthArr: number[] = get(db.depths, [sample, segment], []);
    const meanCov = meanCoverage(depthArr, 1, depthArr.length).toFixed(2);
    const medianCov = medianCoverage(depthArr, 1, depthArr.length).toFixed(2);
    const genomeCov = genomeCoverage(depthArr, 1, depthArr.length, db.chartOptions.low_coverage_threshold).toFixed(2);
    let coverageDepth;
    if (position <= depthArr.length) {
      coverageDepth = depthArr[position - 1].toLocaleString();
    } else {
      coverageDepth = "N/A";
    }
    const row = [
      sample,
      coverageDepth,
      segment,
      depthArr.length > 0 ? `1-${depthArr.length}` : `${depthArr.length}-${depthArr.length}`,
      depthArr.length,
      meanCov,
      medianCov,
      genomeCov
    ];
    rows.push(row);
  });
  return rows;
}

export function tooltipFormatter(db: WgsCovPlotDB) {
  const selectedSamples = db.chartOptions.selectedSamples;
  if (isNil(db.segments)) {
    const depths: SampleDepths = unwrap(db.depths) as SampleDepths;
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      const output = "";
      const [{
        axisIndex,
        axisValue: position
      }] = params;
      if (axisIndex >= selectedSamples.length || db.variants === undefined) {
        return output;
      }
      const sample = selectedSamples[axisIndex];
      const sampleDepths = depths[sample];
      // change tooltip position if it's a click event and the tooltip is not already showing
      if (!db.tooltipOptions.show) {
        setState("tooltipOptions", "left", db.tooltipOptions.x + 10);
        setState("tooltipOptions", "top", db.tooltipOptions.y);
      }
      // set tooltip state
      setState("tooltipOptions", "show", true);
      setState("tooltipOptions", "sample", sample);
      setState("tooltipOptions", "position", position);
      setState("tooltipOptions", "depth", sampleDepths[position - 1]);

      const [dataZoom] = db.chart.getOption().dataZoom;
      let zoomStart = 1
      let zoomEnd = db.ref_seq?.length;
      if (dataZoom !== undefined) {
        zoomStart = Math.floor(dataZoom.startValue);
        zoomEnd = Math.floor(dataZoom.endValue);
      }
      if (isNil(zoomEnd)) {
        return;
      }
      let positionRows: (string | number)[][] = [];
      const tables: { headers: string[], rows: (string | number)[][] }[] = [];
      const samples = (db.chartOptions.crossSampleComparisonInTooltips) ? db.chartOptions.selectedSamples : [sample];
      positionRows = getVariantComparison(db, position, samples);
      tables.push({headers: ["Position Info", ""], rows: positionRows});
      if (!isEmpty(positionRows)) {
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getCoverageStatComparison(db, zoomStart, zoomEnd, position);
          } else {
            const meanCov = meanCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            const medianCov = medianCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            const genomeCov = genomeCoverage(sampleDepths, zoomStart, zoomEnd, db.chartOptions.low_coverage_threshold).toFixed(2);
            coverageStatRows = [
              [
                "Range",
                `${zoomStart.toLocaleString()} - ${zoomEnd.toLocaleString()}`,
              ],
              ["Mean Coverage", `${meanCov}X`],
              ["Median Coverage", `${medianCov}X`],
              [`Genome Coverage (>= ${db.chartOptions.low_coverage_threshold}X)`, `${genomeCov}%`],
            ];
          }
          tables.push({
            headers: ["Coverage View Stats", ""],
            rows: coverageStatRows,
          });
        }
      }
      setState("tooltipOptions", "tables", tables);
      return;
    };
  } else {
    const depths: SampleSegmentDepths = unwrap(db.depths) as SampleSegmentDepths;
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      const [{axisIndex}] = params;
      let [{axisValue: position}] = params;
      let positionRows: (string | number)[][] = [];
      const tables = [];
      const sample: string = selectedSamples[axisIndex];
      const segment: string = whichSegment(position, db);
      const sequence: string = get(db, ["segments_ref_seq", sample, segment], "");
      const segmentLength = sequence.length;
      // convert to pos in segment
      const start = get(db, ["segCoords", segment, "start"], Number.MAX_VALUE);
      position = position - start + 1;
      // eslint-disable-next-line
      let coverageDepth: any;
      const sampleSegDepths: number[] = get(depths, [sample, segment]);
      const refID: string = get(db, ["segments_ref_id", sample, segment], "");
      if (segmentLength === 0) {
        coverageDepth = "N/A";
      } else {
        if (position <= segmentLength) {
          // get coverage depth for pos
          coverageDepth = sampleSegDepths[position - 1].toLocaleString();
        } else {
          coverageDepth = `No sequence at this position. 
                        Reference sequence ${refID} 
                        is only ${sampleSegDepths.length} bp`;
        }
      }
      if (!db.tooltipOptions.show) {
        setState("tooltipOptions", "left", db.tooltipOptions.x + 10);
        setState("tooltipOptions", "top", db.tooltipOptions.y);
      }
      // set tooltip state
      setState("tooltipOptions", "show", true);
      setState("tooltipOptions", "sample", sample);
      setState("tooltipOptions", "position", position);
      setState("tooltipOptions", "depth", coverageDepth);
      const samples = (db.chartOptions.crossSampleComparisonInTooltips) ? db.chartOptions.selectedSamples : [sample];
      positionRows = getSegmentVariantComparison(db, samples, segment, position)
      tables.push({headers: ["Position Info", ""], rows: positionRows});
      if (positionRows.length) { // write rows to table
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getSegmentCoverageStatComparison(db, segment, position);
          } else {
            const meanCov = meanCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            const medianCov = medianCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            const genomeCov = genomeCoverage(sampleSegDepths, 1, segmentLength, db.chartOptions.low_coverage_threshold).toFixed(2);
            coverageStatRows = [
              ["Mean Coverage", `${meanCov}X`],
              ["Median Coverage", `${medianCov}X`],
              [`Genome Coverage (>= ${db.chartOptions.low_coverage_threshold}X)`, `${genomeCov}%`],
            ];
          }
          tables.push({headers: ["Coverage View Stats", ""], rows: coverageStatRows});
        }
        if (!isEmpty(db.primer_matches)) {
          const primerInfoRows = getPrimerInfo(sample, position, segment, db)
          if (primerInfoRows.length) {
            tables.push({headers: ["Primer Info", ""], rows: primerInfoRows});
          }
        }
      }
      setState("tooltipOptions", "tables", tables);
      return;
    }
  }
}

export const getTooltips = (db: WgsCovPlotDB) => {
  return [
    {
      trigger: "axis",
      enterable: true,
      triggerOn: db.tooltipOptions.showTooltip ? db.chartOptions.tooltipTriggerOn : "none",
      appendToBody: false,
      renderMode: "html",
      showContent: true,
      confine: true,
      position: "cursor",
      className: "hidden",
      axisPointer: {
        type: "line"
      },
      formatter: tooltipFormatter(db),
    },
  ];
}