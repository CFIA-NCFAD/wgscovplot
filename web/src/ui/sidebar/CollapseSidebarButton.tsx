import {Component} from "solid-js";
import {setState, state} from "../../state";

export const CollapseSidebarButton: Component = () => {
  return <div class="absolute right-2 top-2 cursor-pointer dark:text-gray-300"
              title={state.chartOptions.sidebarCollapsed ? "Click to expand sidebar" : "Click to collapse sidebar"}
              onClick={() => setState("chartOptions", "sidebarCollapsed", !state.chartOptions.sidebarCollapsed)}>
    <svg xmlns="http://www.w3.org/2000/svg" class="h-6 w-6" fill="none" viewBox="0 0 24 24"
         stroke="currentColor">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
            d={state.chartOptions.sidebarCollapsed ? "M13 5l7 7-7 7M5 5l7 7-7 7" : "M11 19l-7-7 7-7m8 14l-7-7 7-7"}/>
    </svg>
  </div>
}