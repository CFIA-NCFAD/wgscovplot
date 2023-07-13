import {Component} from "solid-js";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const CovPlotYScale: Component = () => {
  return <div class="mt-2">
    <fieldset class="">
      <legend class="mb-1">Y-axis scale</legend>
      <input id="y-scale-log"
             class="peer/log form-radio mr-2 border-slate-700 "
             type="radio"
             name="log"
             checked={state.chartOptions.scaleType === "log"}
             onChange={() => setState("chartOptions", "scaleType", "log")}/>
      <label for="y-scale-log" class="peer-checked/log:text-sky-500">Log</label>
      <HelpIcon
        helpMsg="Logarithmic scale is recommended for coverage plots especially if you want to show lower coverage regions alongside high coverage. Otherwise, the high coverage regions drown out the low coverage regions."/>
      <input id="y-scale-linear"
             class="peer/linear form-radio ml-4 mr-2 border-slate-700"
             type="radio"
             name="linear"
             checked={state.chartOptions.scaleType === "value"}
             onChange={() => setState("chartOptions", "scaleType", "value")}/>
      <label for="y-scale-linear" class="peer-checked/linear:text-sky-500">Linear</label>
    </fieldset>

  </div>
}