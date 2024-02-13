import {Component, Show} from "solid-js";
import {state} from "../../state";
import {ToggleLowCoverageHighlight} from "./ToggleLowCoverageHighlight";
import {LowCovColour} from "./LowCovColour";
import {LowCovThreshold} from "./LowCovThreshold";
import {LowCovThresholdLine} from "./LowCovThresholdLine";

export const LowCovRegionsOptions: Component = () => {
  return <div class="border border-gray-300 rounded dark:border-slate-700 p-1 pl-2 pb-1 mt-1 -ml-1">
    <p class="text-gray-700 font-bold text-lg dark:text-gray-400">
      {`Low coverage regions (< ${state.chartOptions.low_coverage_threshold}X)`}
    </p>
    <ToggleLowCoverageHighlight/>
    <Show when={state.chartOptions.showLowCovRegions}>
      <LowCovColour/>
      <LowCovThreshold/>
      <LowCovThresholdLine/>
    </Show>
  </div>
}