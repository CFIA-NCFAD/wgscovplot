import {Component, For} from "solid-js";
import {HelpIcon} from "../components/HelpIcon";
import {setState, state} from "../../state";

export const NuclColourSelect: Component = () => {
  return <div class="mt-2">
    <p class="dark:text-gray-300">
      Variant colour by nucleotide
      <HelpIcon
        helpMsg="Select a colour for each nucleotide when displaying variant call results overlaid on the coverage subplots."/>
    </p>
    <For each={Object.keys(state.chartOptions.ntColor)}>
      {(nt: string) => {
        return <div class="inline-block">
          <label class="mr-2 dark:text-gray-300" for={nt}>{nt}</label>
          <input type="color" id={nt} class="hover:ring h-5 w-5 mr-2"
                 value={state.chartOptions.ntColor[nt] as string}
                 onChange={(e) => setState("chartOptions", "ntColor", nt, e.currentTarget.value)}/>
        </div>
      }}
    </For>
  </div>
}