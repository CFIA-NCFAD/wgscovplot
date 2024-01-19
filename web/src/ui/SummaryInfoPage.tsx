import {Component} from "solid-js";
import {state} from "../state";
import {isNil, get} from "lodash";
import {MosdepthInfo, SampleMosdepthInfo, SampleSegmentMosdepthInfo} from "../db";
import {TooltipTable} from "./components/TooltipTable";

const COLUMNS = [
  ["n_zero_coverage", "# 0X Coverage Positions", "integer"],
  ["zero_coverage_coords", "0X Coverage Regions", "string"],
  ["low_coverage_threshold", "Low Coverage Threshold (X)", "integer"],
  ["n_low_coverage", "Low Coverage Positions", "integer"],
  ["low_coverage_coords", "Low Coverage Regions", "string"],
  ["genome_coverage", "Genome Coverage (%)", "percent"],
  ["mean_coverage", "Mean Coverage Depth (X)", "float"],
  ["median_coverage", "Median Coverage Depth (X)", "integer"],
  ["ref_seq_length", "Ref Sequence Length (bp)", "integer"],
  ["max_depth", "Max Depth (X)", "integer"],
]


function buildRow(depthInfo: MosdepthInfo, row: string[]) {
  for (const [field, , dtype] of COLUMNS) {
    let val = get(depthInfo, field, "");
    if (dtype === "integer") {
      val = parseInt(val).toLocaleString();
    } else if (dtype === "float") {
      val = parseFloat(val).toFixed(2);
    } else if (dtype === "percent") {
      val = (parseFloat(val) * 100).toFixed(2) + "%";
    }
    row.push(val);
  }
}

export const SummaryInfoPage: Component = () => {
  let headers: string[];
  const rows: string[][] = [];

  if (isNil(state.segments)) {
    const sampleMosdepthInfos = state.mosdepth_info as SampleMosdepthInfo;
    if (isNil(sampleMosdepthInfos)) {
      return null;
    }
    headers = ["Sample", ...COLUMNS.map(([, y]) => y)];
    for (const [sample, depthInfo] of Object.entries(sampleMosdepthInfos)) {
      const row: string[] = [sample];
      buildRow(depthInfo, row);
      rows.push(row);
    }
  } else {
    const sampleSegDepthInfos = state.mosdepth_info as SampleSegmentMosdepthInfo;
    if (isNil(sampleSegDepthInfos)) {
      return null;
    }
    headers = [
      "Sample",
      "Segment",
      ...COLUMNS.map(([, y]) => y)
    ];
    for (const [sample, segDepthInfo] of Object.entries(sampleSegDepthInfos)) {
      for (const [segment, depthInfo] of Object.entries(segDepthInfo)) {
        const row: string[] = [sample, segment];
        buildRow(depthInfo, row);
        rows.push(row);
      }
    }
  }
  return (
      <div class="w-full">
        <TooltipTable headers={headers} rows={rows}/>
      </div>
    );
}