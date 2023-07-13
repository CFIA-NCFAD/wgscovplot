import {Component} from "solid-js";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const LowCovThresholdLine: Component = () => {
  return <div class="w-full mt-2">
    <label for="low-cov-threshold-line-width" class="mr-2">Line width</label>
    <input type="number" min="0" step="0.05" id="low-cov-threshold-line-width"
           value={state.chartOptions.lowCovThresholdLineWidth}
           class="border border-gray-300 rounded justify-text-right px-1 w-1/6"
           onChange={(e) => setState("chartOptions", "lowCovThresholdLineWidth", parseFloat(e.currentTarget.value))}></input>
    <div class="m-2"></div>
    <label for="low-cov-threshold-line-colour" class="mt-2">Line colour</label>
    <HelpIcon helpMsg="Set the colour of low coverage threshold line in the coverage plots."/>
    <input type="color" value={state.chartOptions.lowCovThresholdLineColour}
           class="hover:ring h-5 w-5" id="low-cov-threshold-line-colour"
           onChange={(e) => setState("chartOptions", "lowCovThresholdLineColour", e.currentTarget.value)}/>
  </div>
}