import {Component, Show} from "solid-js";
import {CollapsiblePanel} from "./CollapsiblePanel";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";
import {TooltipTriggerSelect} from "./TooltipTriggerSelect";

export const TooltipOptionsPanel: Component = () => {
  return (
    <CollapsiblePanel title="Tooltip Options" defaultCollapsed={true}>
      <div class="mt-2">
        <input type="checkbox" id="show-tooltip" class="hover:ring h-4 w-4 form-check"
               checked={state.tooltipOptions.showTooltip}
               onChange={(e) => setState("tooltipOptions", "showTooltip", e.currentTarget.checked)}/>
        <label class="ml-1" for="show-tooltip">Show tooltip?</label>
        <HelpIcon
          helpMsg="Show a tooltip on mouse over or click with information about a position/variant in a sample."/>
      </div>
      <Show when={state.tooltipOptions.showTooltip}>
        <TooltipTriggerSelect/>
        <div class="mt-2">
          <input type="checkbox" id="show-cross-sample-comp-tooltips" class="hover:ring h-4 w-4 form-check"
                 checked={state.chartOptions.crossSampleComparisonInTooltips}
                 onChange={(e) => setState("chartOptions", "crossSampleComparisonInTooltips", e.currentTarget.checked)}/>
          <label class="ml-1" for="show-cross-sample-comp-tooltips">Show sample comparison info in tooltips?</label>
          <HelpIcon helpMsg="Show a comparison of the selected sample to the other samples in the tooltip."/>
        </div>
        <div class="mt-2">
          <input type="checkbox" id="show-tooltips-variant-sites-only" class="hover:ring h-4 w-4 form-check"
                 checked={state.tooltipOptions.variantSitesOnly}
                 onChange={(e) => setState("tooltipOptions", "variantSitesOnly", e.currentTarget.checked)}/>
          <label class="ml-1" for="show-tooltips-variant-sites-only">Show tooltips only for variant sites?</label>
          <HelpIcon helpMsg="Show a tooltip only for variant sites (i.e. sites with a non-reference allele)."/>
        </div>
        <div class="mt-2">
          <label for="tooltip-top-position">Tooltip top position</label>
          <HelpIcon helpMsg="Set the top position of the tooltip."/>
          <input id="tooltip-top-position" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
                 type="number" min="0" step="1" value={state.tooltipOptions.top}
                 onChange={(e) => setState("tooltipOptions", "top", parseFloat(e.currentTarget.value))}/>
        </div>
        <div class="mt-2">
          <label for="tooltip-left-position">Tooltip left position</label>
          <HelpIcon helpMsg="Set the left position of the tooltip."/>
          <input id="tooltip-left-position" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
                 type="number" min="0" step="1" value={state.tooltipOptions.left}
                 onChange={(e) => setState("tooltipOptions", "left", parseFloat(e.currentTarget.value))}/>
        </div>
        <div class="mt-2">
          <button class="btn rounded
          dark:bg-slate-900
          dark:hover:bg-slate-700
          hover:bg-gray-200
          bg-slate-300
          text-gray-900
          dark:text-gray-100
          p-1"
            onClick={() => {
            // get height and width of screen
            const height = window.innerHeight;
            const width = window.innerWidth;
            // set tooltip position to top middle of screen
            setState("tooltipOptions", "top", height * 0.1);
            setState("tooltipOptions", "left", width * 0.5);
          }}>
            Reset tooltip position
          </button>
        </div>
      </Show>
      <code class="text-xs">
        x,y: {state.tooltipOptions.x},{state.tooltipOptions.y}
        <br/>
        top,left: {state.tooltipOptions.top},{state.tooltipOptions.left}
      </code>
    </CollapsiblePanel>
  );
}