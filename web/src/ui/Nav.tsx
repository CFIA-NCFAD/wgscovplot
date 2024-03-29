import {Component, Show} from "solid-js";
import {setState, state} from "../state";
import {ThemeButton} from "./ThemeButton";
import {isNil} from "lodash";

export const Nav: Component = () => {
  return <nav
    class="relative py-1 bg-gray-100 dark:bg-gray-800 text-gray-500 dark:text-gray-300 shadow navbar">
    <div class="w-full px-4 flex flex-wrap">
      <div class="items-center justify-between">
        <a
          class="font-bold leading-relaxed inline-block mr-4 py-1 whitespace-nowrap uppercase text-gray-700 dark:text-gray-300"
          href="#">wgscovplot</a>
        <a
          class={"text-sm inline-block mr-4 py-1 " + (state.activePage === "chart" ? " font-bold text-gray-700 dark:text-gray-300" : "text-gray-500 dark:text-gray-400")}
          href="#" onClick={() => setState("activePage", "chart")}>
          Coverage Plot
        </a>
        <Show when={isNil(state.segments)}>
          <a
            class={"text-sm inline-block mr-4 py-1 " + (state.activePage === "heatmap" ? " font-bold text-gray-700 dark:text-gray-300" : "text-gray-500 dark:text-gray-400")}
            href="#" onClick={() => setState("activePage", "heatmap")}>
            Variant Heatmap
          </a>
        </Show>
        <a
          class={"text-sm inline-block mr-4 py-1 " + (state.activePage === "info" ? " font-bold text-gray-700 dark:text-gray-300" : "text-gray-500 dark:text-gray-400")}
          href="#" onClick={() => setState("activePage", "info")}>
          Summary Information
        </a>
        <a
          class={"text-sm inline-block mr-4 py-1 " + (state.activePage === "about" ? " font-bold text-gray-700 dark:text-gray-300" : "text-gray-500 dark:text-gray-400")}
          href="#" onClick={() => setState("activePage", "about")}>
          About
        </a>
      </div>
      <div class="flex-grow"></div>
      <div class="flex-shrink">
        <ThemeButton/>
      </div>
    </div>
  </nav>
}