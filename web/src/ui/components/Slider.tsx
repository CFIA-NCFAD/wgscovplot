import {Show} from "solid-js";
import {isNil} from "lodash";

import type {Component} from "solid-js";

import {HelpIcon} from "./HelpIcon";

type SliderProps = {
  label: string,
  min: number,
  max: number,
  step: number,
  value: number,
  helpInfo?: string,
  onChange: (e: never) => void
}
export const Slider: Component<SliderProps> = (props) => {
  return <div>
    <div class="w-full mt-2 inline-block">
      <label>{props.label}</label>
      <Show when={!isNil(props.helpInfo)}>
        <HelpIcon helpMsg={isNil(props.helpInfo) ? props.label : props.helpInfo}/>
      </Show>
      <span class="absolute right-4 text-gray-500 text-sm pr-1">
        {props.value}%
      </span>
    </div>
    <input type="range" min={props.min} max={props.max} step={props.step}
           value={props.value}
           class="w-full"
           // @ts-expect-error ignore TS
           onChange={props.onChange}></input>
  </div>
}
