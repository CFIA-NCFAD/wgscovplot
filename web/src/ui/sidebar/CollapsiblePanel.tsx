import {Component, createSignal} from "solid-js";

export const CollapsiblePanel: Component<{ title: string, defaultCollapsed?: boolean, children: never }> = (props) => {
  const [isOpen, setIsOpen] = createSignal(!props.defaultCollapsed);
  return <div class="border border-gray-300 rounded dark:border-slate-700 w-full dark:text-gray-300">
    <button class="hover:bg-gray-300 dark:hover:bg-slate-900 flex w-full justify-between " onClick={() => setIsOpen(!isOpen())}>
      <p class="justify-self-start text-gray-700 font-bold text-lg dark:text-gray-400 flex-none ml-1">
        {props.title}
      </p>
      <svg xmlns="http://www.w3.org/2000/svg" class="h-6 w-6 flex-none pr-2" fill="none" viewBox="0 0 24 24"
           stroke="currentColor">


        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="3"
              d={isOpen() ? "M 13 11 L 23 21 M 13 11 L 3 21" : "M 13 21 L 3 11 M 13 21 L 23 11"}/>
      </svg>
      <span class="hidden justify-self-end flex-none">{isOpen() ? "Hide" : "Show"}</span>
    </button>
    <div class={"p-1 pl-2 " + (isOpen() ? "" : "hidden")}>
      {props.children}
    </div>
  </div>
}