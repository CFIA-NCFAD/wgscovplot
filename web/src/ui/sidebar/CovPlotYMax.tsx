import {Component} from "solid-js";
import {HelpIcon} from "../components/HelpIcon";
import {setState, state} from "../../state";

export const CovPlotYMax: Component = () => {
  return <div>
    <div class="w-full mt-2 inline-block">
      <label for="y-max">Y-max</label>
      <HelpIcon
        helpMsg="Adjust the Y-axis coverage depth max or limit value. By default, it's dependent on samples being shown."/>
      <input type="number" min="0" step="1" id="y-max"
             value={state.chartOptions.yMax}
             class="ml-3 border border-gray-300 rounded justify-text-right px-1 w-1/4"
             onChange={(e) => setState("chartOptions", "yMax", parseInt(e.currentTarget.value))}></input>
    </div>
  </div>
}