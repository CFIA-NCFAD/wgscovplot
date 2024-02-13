import {Component, For} from "solid-js";
import {createDraggable, transformStyle} from "@thisbeyond/solid-dnd";
import {setState, state} from "../state";
import {HelpIcon} from "./components/HelpIcon";
import {TooltipTable} from "./components/TooltipTable";

export const ChartTooltip: Component = () => {
  const draggable = createDraggable(1);

  // get the index of the column that contains the sample name
  function getHighlightColumnIndex(rows: (string | number)[][], sample: string) {
    for (let i = 0; i < rows.length; i++) {
      for (let j = 1; j < rows[i].length; j++) {
        if (rows[i][j] === sample) {
          return j;
        }
      }
    }
    return undefined;
  }

  return <div ref={draggable.ref}
              class={"absolute bg-gray-100 dark:bg-gray-800 shadow border border-gray-300 dark:border-gray-800 rounded p-1 min-w-[10%] " + (state.tooltipOptions.show ? "" : "hidden")}
              style={{
                ...{top: state.tooltipOptions.top + "px", left: state.tooltipOptions.left + "px"},
                ...transformStyle(draggable.transform)
              }}>
    <div class="w-full text-sm text-gray-800 dark:text-gray-300 dark:bg-gray-800">
      <button class="absolute right-1 hover:bg-slate-300 py-1 px-2 rounded text-2xl"
              onClick={() => setState("tooltipOptions", "show", false)}>
        ⊗
      </button>

      <div class="text-lg font-bold " style={{cursor: "move"}} {...draggable.dragActivators}>
        {state.tooltipOptions.sample}
        <HelpIcon helpMsg="Drag the tooltip to move it around."/>
      </div>

      <p class="ml-2">Position: {state.tooltipOptions.position.toLocaleString()}</p>
      <p class="ml-2">Depth: {state.tooltipOptions.depth.toLocaleString()}</p>
      <For each={state.tooltipOptions.tables}>
        {(table) => {
          const props: { highlightColumn?: number } = {};
          props.highlightColumn = getHighlightColumnIndex(table.rows, state.tooltipOptions.sample);
          return <div class="mt-4">
            <TooltipTable headers={table.headers}
                          rows={table.rows}
                          {...props}
            />
          </div>;
        }}
      </For>
    </div>
  </div>
}