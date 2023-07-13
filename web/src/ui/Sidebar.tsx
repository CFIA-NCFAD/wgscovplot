import {state} from "../state";
import type {Component} from 'solid-js';
import {Show} from "solid-js";
import {isNil} from "lodash";
import {SampleSelect} from "./sidebar/SampleSelect";
import {TooltipOptionsPanel} from "./sidebar/TooltipOptionsPanel";
import {AxisOptionsPanel} from "./sidebar/AxisOptionsPanel";
import {MutationOptionsPanel} from "./sidebar/MutationOptionsPanel";
import {DisplayOptionsPanel} from "./sidebar/DisplayOptionsPanel";
import {SegmentSelect} from "./sidebar/SegmentSelect";
import {CollapseSidebarButton} from "./sidebar/CollapseSidebarButton";


export const Sidebar: Component = () => {
  return <aside
    class={"bg-gray-100 dark:bg-slate-800 overflow-auto " + (!state.chartOptions.sidebarCollapsed ? "sm:w-1/5 " : "")}>
    <div class="sticky top-0 pt-4 px-4 w-full text-sm">
      <CollapseSidebarButton/>
      <div class={state.chartOptions.sidebarCollapsed ? "hidden" : ""}>
        <SampleSelect/>
        <Show when={!isNil(state.segments)}>
          <SegmentSelect/>
        </Show>
        <AxisOptionsPanel/>
        <Show when={!isNil(state.variants)}>
          <MutationOptionsPanel/>
        </Show>
        <TooltipOptionsPanel/>
        <DisplayOptionsPanel/>
      </div>
    </div>
  </aside>
}

