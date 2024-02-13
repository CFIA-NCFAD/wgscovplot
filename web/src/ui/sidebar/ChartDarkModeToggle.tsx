// dark mode toggle
import {Component} from "solid-js";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";

export const ChartDarkModeToggle: Component = () => {
  return <div class="mt-2">
    <input type="checkbox"
           id="dark-mode-toggle"
           class="hover:ring h-4 w-4 form-check"
           checked={state.chartOptions.darkMode}
           onChange={(e) => setState("chartOptions", "darkMode", e.currentTarget.checked)}
    />
    <label class="ml-1" for="dark-mode-toggle">Chart Dark Mode?</label>
    <HelpIcon helpMsg="Toggle ECharts dark mode."/>
  </div>
}