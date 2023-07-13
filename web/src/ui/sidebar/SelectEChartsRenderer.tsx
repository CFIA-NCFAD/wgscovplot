// radio input for selecting ECharts renderer
import {Component} from "solid-js";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const SelectEChartsRenderer: Component = () => {
  return <div class="mt-2">
    <fieldset class="">
      <legend class="mb-1">ECharts renderer</legend>
      <input id="echarts-renderer-canvas"
             class="peer/canvas form-radio mr-2 border-slate-700 "
             type="radio"
             name="renderer"
             checked={state.chartOptions.renderer === "canvas"}
             onChange={() => setState("chartOptions", "renderer", "canvas")}/>
      <label for="echarts-renderer-canvas" class="peer-checked/canvas:text-sky-500">Canvas</label>
      <HelpIcon
        helpMsg="Canvas renderer is recommended for large coverage plots and when high performance is required. Outputs file in PNG format."/>
      <input id="echarts-renderer-svg"
             class="peer/svg form-radio ml-4 mr-2 border-slate-700"
             type="radio"
             name="renderer"
             checked={state.chartOptions.renderer === "svg"}
             onChange={() => setState("chartOptions", "renderer", "svg")}/>
      <label for="echarts-renderer-svg" class="peer-checked/svg:text-sky-500">SVG</label>
      <HelpIcon
        helpMsg="SVG renderer is recommended for publication quality graphics. Outputs file in scale vector graphic (SVG) format. Performance may be slightly poorer than Canvas rendering."/>
    </fieldset>
  </div>
}