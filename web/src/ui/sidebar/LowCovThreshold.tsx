import {Component} from "solid-js";
import {setState, state} from "../../state";

export const LowCovThreshold: Component = () => {
  return <div class="w-full mt-2 inline-block">
    <label for="low-cov-threshold" class="mr-2">Threshold</label>
    <input type="number" min="0" step="1" id="low-cov-threshold"
           value={state.chartOptions.low_coverage_threshold}
           class="border border-gray-300 rounded justify-text-right px-1 w-1/6"
           onChange={(e) => setState("chartOptions", "low_coverage_threshold", parseInt(e.currentTarget.value))}></input>
    <span class="ml-1">X</span>
  </div>
}