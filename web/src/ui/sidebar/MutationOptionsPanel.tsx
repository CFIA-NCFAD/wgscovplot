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
      </Show>
    </Show>
  </CollapsiblePanel>
}