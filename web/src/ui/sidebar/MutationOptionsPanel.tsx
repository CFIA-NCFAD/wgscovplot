import {Component, Show} from "solid-js";
import {CollapsiblePanel} from "./CollapsiblePanel";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const MutationOptionsPanel: Component = () => {
  return <CollapsiblePanel title="Mutation Options" defaultCollapsed={true}>
    <div class="mt-2">
      <input type="checkbox" id="show-variants" class="hover:ring h-4 w-4 form-check"
             checked={state.chartOptions.showVariants}
             onChange={(e) => setState("chartOptions", "showVariants", e.currentTarget.checked)}/>
      <label class="ml-1" for="show-variants">Show variants?</label>
      <HelpIcon
        helpMsg="Show variant calling/mutation information in the subplots. Variant calls need to be shown in order for coverage plot tooltips with variant calling results to be shown on click or mouse over."/>
    </div>
    <Show when={state.chartOptions.showVariants} keyed>
      <div class="mt-2">
        <label class="ml-1" for="variant-bar-width">Variant bar width</label>
        <HelpIcon helpMsg="Set width of variant bars in the subplots."/>
        <input type="number" id="variant-bar-width" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
               value={state.chartOptions.variantBarWidth} step="0.1"
               onChange={(e) => setState("chartOptions", "variantBarWidth", parseFloat(e.currentTarget.value))}/>
      </div>
      <div class="mt-2">
        <input type="checkbox" id="show-mutation-labels" class="hover:ring h-4 w-4 form-check"
               checked={state.chartOptions.showVariantLabels}
               onChange={(e) => setState("chartOptions", "showVariantLabels", e.currentTarget.checked)}/>
        <label class="ml-1" for="show-mutation-labels">Show mutation/variant calling labels in x-axis?</label>
        <HelpIcon helpMsg="Show mutation labels in the subplots."/>
      </div>
      <Show when={state.chartOptions.showVariantLabels} keyed>
        <div class="mt-2">
          <input type="checkbox" id="hide-overlapping-variant-labels" class="hover:ring h-4 w-4 form-check"
                 checked={state.chartOptions.hideOverlappingVariantLabels}
                 onChange={(e) => setState("chartOptions", "hideOverlappingVariantLabels", e.currentTarget.checked)}/>
          <label class="ml-1" for="hide-overlapping-variant-labels">Hide overlapping variant labels?</label>
          <HelpIcon helpMsg="Hide mutation/variant labels that overlap with each other."/>
        </div>
        <div class="mt-2">
          <label for="subplot-title-font-size">Variant labels rotation angle</label>
          <HelpIcon helpMsg="Set rotation angle for variant labels"/>
          <input id="subplot-title-font-size" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
                 type="number" min="-180" step="1" value={state.chartOptions.variantLabelsRotation}
                 onChange={(e) => setState("chartOptions", "variantLabelsRotation", parseFloat(e.currentTarget.value))}/>
        </div>
      </Show>
    </Show>
  </CollapsiblePanel>
}