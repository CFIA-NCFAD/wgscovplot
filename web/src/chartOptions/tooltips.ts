import {get, isNil, find, uniqBy, flatMap, keys, isEmpty} from "lodash";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {SampleDepths, SegmentCoords, VariantCall, WgsCovPlotDB} from "../db";
import {unwrap} from "solid-js/store";
import {setState} from "../state";
import {whichSegment} from "./segmented/getSegmentsInfo";
import {getPrimerInfo} from "./segmented/pcr";


export const getVariantComparison = (db: WgsCovPlotDB, position: number): Array<Array<string>> => {
  if (isEmpty(db.variants)) {
    return [];
  }
  let rows: any[][] = [];
  let variantArr: VariantCall[] = [];
  for (let sample of db.chartOptions.selectedSamples) {
    // @ts-ignore
    let variant = find(db.variants[sample], {POS:  position}, 0) as VariantCall;
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
    let variantCall = find(db.variants[sample][segment], {POS: position}, 0);
    if (!isNil(variantCall)) {
      variantArr.push(variantCall);
    } else {
      variantArr.push({
        "Sample": sample,
        "POS": position,
        "REF_ID": isNil(db.segments_ref_id[sample][segment]) ? "" : db.segments_ref_id[sample][segment],
        "Segment": segment,
        // @ts-ignore
        "Segment Length": isNil(db.depths[sample][segment]) ? 0 : db.depths[sample][segment].length,
        "REF_SEQ": isNil(db.segments_ref_seq[sample][segment]) ? "" : db.segments_ref_seq[sample][segment][position - 1],
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
    let depthArr: number[] = isNil(db.depths[sample][segment]) ? [] : db.depths[sample][segment];
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
      depthArr.length > 0 ?`1-${depthArr.length}` : `${depthArr.length}-${depthArr.length}`,
      depthArr.length,
      meanCov,
      medianCov,
      genomeCov
    ];
    rows.push(row);
  });
  return rows;
}

export const tooltipFormatter = (db: WgsCovPlotDB) => {
  if (isNil(db.segments)) {
    let depths: SampleDepths = unwrap(db.depths) as SampleDepths;
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      let selectedSamples = db.chartOptions.selectedSamples;
      let output = "";
      let [{
        axisIndex,
        axisValue: position
      }] = params;
      if (axisIndex >= selectedSamples.length || db.variants === undefined) {
        return output;
      }
      let sample = selectedSamples[axisIndex];
      let sampleDepths = depths[sample];
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

      let [dataZoom] = db.chart.getOption().dataZoom;
      let zoomStart = 1
      let zoomEnd = db.ref_seq.length;
      if (dataZoom !== undefined) {
        zoomStart = Math.floor(dataZoom.startValue);
        zoomEnd = Math.floor(dataZoom.endValue);
      }
      let positionRows: string[][] = [];
      let tables = [];
      const isVariantBar = params.find(x => x.componentSubType === "bar");
      if (isVariantBar) {
        if (db.chartOptions.crossSampleComparisonInTooltips) {
          positionRows = getVariantComparison(db, position);
        } else {
          // @ts-ignore
          let foundObj: any = find(db.variants[sample], {POS: position}, 0);
          if (!isNil(foundObj)) {
            for (const [key, value] of Object.entries(foundObj)) {
              positionRows.push([key, value as string]);
            }
          }
        }
        tables.push({headers: ["Variant Info", ""], rows: positionRows});
      } else {
        positionRows = [
          ["Sample", sample],
          ["Position", position.toString()],
          ["Sequence", db.ref_seq[position - 1]]
        ]
        tables.push({headers: ["Position Info", ""], rows: positionRows});
      }
      if (positionRows.length) {
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getCoverageStatComparison(db, zoomStart, zoomEnd, position);
          } else {
            let meanCov = meanCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            let medianCov = medianCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            let genomeCov = genomeCoverage(sampleDepths, zoomStart, zoomEnd, db.chartOptions.low_coverage_threshold).toFixed(2);
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
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      let [{
        axisIndex,
        axisValue: position
      }] = params;
      let positionRows: string[][] = [];
      let tables = [];
      let sample: string = db.chartOptions.selectedSamples[axisIndex];
      let segment: string = whichSegment(position, db);
      let sequence: string = isNil(db.segments_ref_seq[sample][segment]) ? "" : db.segments_ref_seq[sample][segment];
      let segmentLength = sequence.length;
      // convert to pos in segment
      position = position - db.segCoords[segment].start + 1;
      let coverageDepth: any;
      // @ts-ignore
      let sampleSegDepths: number[] = db.depths[sample][segment];
      let refID: string = isNil(db.segments_ref_id[sample][segment]) ? "" : db.segments_ref_id[sample][segment];
      if (segmentLength === 0) {
        coverageDepth = `No result reported for segment ${segment}`;
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
      const isVariantBar = params.find(element => element.componentSubType === "bar");
      if (isVariantBar) { //tooltips at Variant sites
        if (db.chartOptions.crossSampleComparisonInTooltips) {
          positionRows = getSegmentVariantComparison(db, sample, segment, position)
        } else {
          // @ts-ignore
          let foundObj = find(db.variants[sample][segment], {POS: position}, 0);
          if (!isNil(foundObj)) {
            for (const [key, value] of Object.entries(foundObj)) {
              // @ts-ignore
              positionRows.push(...[[key, value]]);
            }
          }
          positionRows.push(["Coverage Depth", coverageDepth]);
        }
        let sortingCriteria = ["Sample", "Segment", "POS", "Segment Length", "Coverage Depth", "REF_ID", "REF_SEQ", "ALT_SEQ", "ALT_FREQ"]
        positionRows.sort((a, b) => sortingCriteria.indexOf(a[0]) - sortingCriteria.indexOf(b[0]))
        tables.push({headers: ["Variant Info", ""], rows: positionRows});
      } else { // tooltips for Non-Variant Sites
        if (position > segmentLength) { //Out of range when segment length < padding array
          positionRows = [
            ["Position", position.toString()],
            [coverageDepth, ""]
          ];
        } // Pos within segment length
        positionRows = [
          ["Sample", sample],
          ["Segment", segment],
          ["POS", position.toString()],
          ["Coverage Depth", coverageDepth],
          ["Segment Length", segmentLength.toString()],
          ["REF_ID", refID],
          ["REF_SEQ", sequence[position - 1]],
          ["ALT_SEQ", ""],
          ["ALT_FREQ", ""],
        ];
        tables.push({headers: ["Position Info", ""], rows: positionRows});
      }
      if (positionRows.length) { // write rows to table
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getSegmentCoverageStatComparison(db, segment, position);
          } else {
            let meanCov = meanCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            let medianCov = medianCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            let genomeCov = genomeCoverage(sampleSegDepths, 1, segmentLength, db.chartOptions.low_coverage_threshold).toFixed(2);
            coverageStatRows = [
              ["Mean Coverage", `${meanCov}X`],
              ["Median Coverage", `${medianCov}X`],
              [`Genome Coverage (>= ${db.chartOptions.low_coverage_threshold}X)`, `${genomeCov}%`],
            ];
          }
          tables.push({headers: ["Coverage View Stats", ""], rows: coverageStatRows});
        }
        if (!isEmpty(db.primer_matches)) {
          let primerInfoRows = getPrimerInfo(sample, position, segment, db)
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