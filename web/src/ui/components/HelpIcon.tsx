import type {Component} from "solid-js";

export const HelpIcon: Component<{ helpMsg: string }> = (props) => {
  return <span
    class="mx-1 border border-blue-400 text-xs text-blue-400 rounded-full cursor-help px-1 font-bold font-mono align-text-top"
    title={props.helpMsg}>?</span>
}