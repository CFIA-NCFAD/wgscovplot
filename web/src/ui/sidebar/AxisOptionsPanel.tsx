import {Component} from "solid-js";
import {CollapsiblePanel} from "./CollapsiblePanel";
import {CovPlotYMax} from "./CovPlotYMax";
import {CovPlotYScale} from "./CovPlotYScale";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const AxisOptionsPanel: Component = () => {
  return <CollapsiblePanel title="Axis Options" defaultCollapsed={true}>
    <CovPlotYMax/>
    <CovPlotYScale/>
    <div class="mt-2">
      <input type="checkbox" id="show-x-axis-labels" class="hover:ring h-4 w-4 form-check"
             checked={state.chartOptions.showXAxisLabel}
             onChange={(e) => setState("chartOptions", "showXAxisLabel", e.currentTarget.checked)}/>
      <label class="ml-1" for="show-x-axis-labels">Show x-axis labels?</label>
      <HelpIcon helpMsg="Show the x-axis labels for genome position in the subplots."/>
    </div>
  </CollapsiblePanel>
}