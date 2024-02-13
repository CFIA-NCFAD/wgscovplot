import {Component} from "solid-js";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const ToggleLowCoverageHighlight: Component = () => {
  return <div class="w-full mt-2 inline-block">
    <div class="form-check form-switch">
      <input type="checkbox"
             checked={state.chartOptions.showLowCovRegions}
             class=" hover:ring h-4 w-4 form-check"
             role="switch"
             id="showLowCovRegions"
             onChange={(e) => setState("chartOptions", "showLowCovRegions", e.currentTarget.checked)}/>
      <label for="showLowCovRegions" class="ml-1">Highlight?</label>
      <HelpIcon helpMsg="Toggle highlighting of low coverage regions in the coverage plots."/>
    </div>
  </div>
}