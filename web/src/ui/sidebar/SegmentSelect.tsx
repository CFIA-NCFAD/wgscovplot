import {Component} from "solid-js";
import {isNil} from "lodash";
import {setState, state} from "../../state";
import {createOptions, Select} from "@thisbeyond/solid-select";
import {HelpIcon} from "../components/HelpIcon";

export const SegmentSelect: Component = () => {
  if (isNil(state.segments)) return null;
  const selectProps = createOptions(state.segments, {
    disable: (value) => state.chartOptions.selectedSegments !== undefined ? state.chartOptions.selectedSegments.includes(value) : false
  });

  return <div class="mb-3">
    <label class="text-gray-700 dark:text-gray-200">
      Segments
    </label>
    <HelpIcon helpMsg="Select genome segments to show in the plot."/>
    <div class="flex-1 flex justify-left items-center gap-5 dark:text-gray-300 dark:bg-slate-600 relative">
      <Select {...selectProps}
              class="bg-white w-full dark:bg-slate-900"
              initialValue={ /*@once*/ state.segments}
              multiple={true}
              onChange={(e) => setState("chartOptions", "selectedSegments", e)}/>
    </div>
  </div>
}