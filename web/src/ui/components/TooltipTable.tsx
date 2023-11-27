import type {Component} from "solid-js";
import {For} from "solid-js";

export const TooltipTable: Component<{ headers: string[], rows: string[][], highlightColumn?: number }> = (props) => {
  return (
    <table class="table-auto p-2 m-1 dark:text-gray-300 dark:bg-slate-800">
      <thead>
      <tr class="border-b border-gray-700">
        <For each={props.headers}>
          {(header, index) => <th class="px-2">{header}</th>}
        </For>
      </tr>
      </thead>
      <tbody>
      <For each={props.rows}>
        {(cells, i) =>
          <tr class={(i() % 2 === 0) ? "bg-gray-300 dark:bg-gray-900" : ""}>
            <For each={cells}>
              {(cell, j) =>
                <td class={"px-2 " +
                  (j() === 0 ? "font-semibold" : "text-right font-mono") +
                  (j() === props.highlightColumn ? " bg-amber-400/25" : "")
                }>
                  {cell}
                </td>
              }
            </For>
          </tr>
        }
      </For>
      </tbody>
    </table>
  );
}
