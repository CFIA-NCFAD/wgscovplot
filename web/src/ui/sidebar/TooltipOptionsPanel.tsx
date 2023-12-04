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
      </Show>
    </CollapsiblePanel>
  );
}