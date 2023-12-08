import {Component, For} from "solid-js";
import {state} from "../state";
import {isNil} from "lodash";

function sortRow(obj: any) {
  let row = Object.keys(obj).map((key) => [key, obj[key]]);
  let sortedRow: any[][];
  if (isNil(state.segments)) {
    let sortingCriteria = ["sample", "n_zero_coverage", "zero_coverage_coords", "n_low_coverage",
      "low_coverage_coords", "genome_coverage", "low_coverage_threshold", "mean_coverage", "median_coverage", "ref_seq_length", "max_depth"]
    sortedRow = row.sort((a, b) => sortingCriteria.indexOf(a[0]) - sortingCriteria.indexOf(b[0]))
  } else {
    let sortingSegmentCriteria = ["sample", "segment", "n_zero_coverage", "zero_coverage_coords", "n_low_coverage",
      "low_coverage_coords", "genome_coverage", "low_coverage_threshold", "mean_coverage", "median_coverage", "ref_seq_length", "max_depth"]
    sortedRow = row.sort((a, b) => sortingSegmentCriteria.indexOf(a[0]) - sortingSegmentCriteria.indexOf(b[0]))
  }
  return sortedRow
}

export const SummaryInfoPage: Component = () => {
  if (isNil(state.segments)) {
    return (
      <div class="w-full">
        <table class="min-w-full text-left text-sm font-light">
          <thead class="border-b bg-white font-medium dark:border-neutral-500 dark:bg-neutral-600">
          <tr>
            <th scope="col" class="px-6 py-4">Sample</th>
            <th scope="col" class="px-6 py-4"># 0 Coverage Positions</th>
            <th scope="col" class="px-6 py-4">0 Coverage Regions</th>
            <th scope="col" class="px-6 py-4"># Low Coverage Positions</th>
            <th scope="col" class="px-6 py-4">Low Coverage Regions</th>
            <th scope="col" class="px-6 py-4">% Genome Coverage</th>
            <th scope="col" class="px-6 py-4">Low Coverage Threshold</th>
            <th scope="col" class="px-6 py-4">Mean Coverage Depth (X)</th>
            <th scope="col" class="px-6 py-4">Median Coverage Depth (X)</th>
            <th scope="col" class="px-6 py-4">Ref Sequence Length (bp)</th>
            <th scope="col" class="px-6 py-4">Max Depth (X)</th>
          </tr>
          </thead>
          <tbody>
          <For each={Object.values(state.mosdepth_info)}>{
            (item, i) =>
              <tr>
                <For each={sortRow(item)}>{
                  (ele, j) =>
                    <td class="px-6 py-4">{ele[1]}</td>
                }</For>
              </tr>
          }</For>
          </tbody>
        </table>
      </div>
    );
  } else {
    return (
      <div class="w-full">
        <table class="min-w-full text-left text-sm font-light">
          <thead class="border-b bg-white font-medium dark:border-neutral-500 dark:bg-neutral-600">
          <tr>
            <th scope="col" class="px-6 py-4">Sample</th>
            <th scope="col" class="px-6 py-4">Segment</th>
            <th scope="col" class="px-6 py-4"># 0 Coverage Positions</th>
            <th scope="col" class="px-6 py-4">0 Coverage Regions</th>
            <th scope="col" class="px-6 py-4"># Low Coverage Positions</th>
            <th scope="col" class="px-6 py-4">Low Coverage Regions</th>
            <th scope="col" class="px-6 py-4">% Genome Coverage</th>
            <th scope="col" class="px-6 py-4">Low Coverage Threshold</th>
            <th scope="col" class="px-6 py-4">Mean Coverage Depth (X)</th>
            <th scope="col" class="px-6 py-4">Median Coverage Depth (X)</th>
            <th scope="col" class="px-6 py-4">Ref Sequence Length (bp)</th>
            <th scope="col" class="px-6 py-4">Max Depth (X)</th>
          </tr>
          </thead>
          <tbody>
          <For each={Object.values(state.mosdepth_info)}>{
            (item, i) =>
              <For each={Object.values(item)}>{
                (item1, j) =>
                  <tr>
                    <For each={sortRow(item1)}>{
                      (ele, k) =>
                        <td class="px-6 py-4">{ele[1]}</td>
                    }</For>
                  </tr>
              }</For>
          }</For>
          </tbody>
        </table>
      </div>
    );
  }
}