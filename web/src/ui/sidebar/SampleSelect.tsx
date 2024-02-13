import {Component} from "solid-js";
import {createOptions, Select} from "@thisbeyond/solid-select";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const SampleSelect: Component = () => {
  const selectProps = createOptions(state.samples, {
    disable: (value) => state.chartOptions.selectedSamples.includes(value)
  });

  return <div class="mb-3 text-gray-700 dark:text-gray-200">
    <label class="">Samples</label>
    <HelpIcon helpMsg="Select samples to show in coverage plot"/>
    <div class="flex-1 flex justify-left items-center gap-5 relative dark:bg-gray-900">
      <Select {...selectProps}
              class="custom w-full bg-white dark:bg-slate-800 dark:text-gray-300 dark:border-gray-900 dark:bg-opacity-0"
              initialValue={ /*@once*/ state.chartOptions.selectedSamples.length === 0 ? state.samples.slice(0, 3) : state.chartOptions.selectedSamples}
              multiple={true}
              onChange={(e) => setState("chartOptions", "selectedSamples", e)}/>
    </div>
  </div>
}