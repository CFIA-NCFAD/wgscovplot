import {Component} from "solid-js";
import {setState, state} from "../../state";

export const TooltipTriggerSelect: Component = () => {
  return <fieldset class="mt-2">
    <legend class="mb-1">Show tooltip on mouse</legend>

    <input id="tooltip-mousemove"
           class="peer/mousemove form-radio mr-2 border-slate-700 "
           type="radio"
           name="mousemove"
           checked={state.chartOptions.tooltipTriggerOn === "mousemove"}
           onChange={() => setState("chartOptions", "tooltipTriggerOn", "mousemove")}/>
    <label for="tooltip-mousemove" class="peer-checked/mousemove:text-sky-500">move</label>

    <input id="tooltip-click"
           class="peer/click form-radio ml-4 mr-2 border-slate-700"
           type="radio"
           name="click"
           checked={state.chartOptions.tooltipTriggerOn === "click"}
           onChange={() => setState("chartOptions", "tooltipTriggerOn", "click")}/>
    <label for="tooltip-click" class="peer-checked/click:text-sky-500">click</label>
  </fieldset>
}