import {Component, For} from "solid-js";
import {CollapsiblePanel} from "./CollapsiblePanel";
import {setState, state} from "../../state";
import {HelpIcon} from "../components/HelpIcon";
import {NuclColourSelect} from "./NuclColourSelect";
import {isNil} from "lodash";

function FeatureColourPickers() {
  return (
    <div class="mt-2">
      <p class="dark:text-gray-300">
        Feature colour
        <HelpIcon
          helpMsg="Select a colour for each feature shown below the coverage subplots."/>
      </p>
      <For each={state.echart_features}>
        {(feature, index) => {
          if (isNil(feature.itemStyle) || isNil(feature.itemStyle.color)) {
            return null;
          }
          if (feature.value.type === "amplicon") {
            return null;
          }
          return (
            <div class="inline-block">
              <label class="mr-2 dark:text-gray-300" for={feature.name}>{feature.name}</label>
              <input type="color" id={feature.name} class="hover:ring h-5 w-5 mr-2"
                     value={feature.itemStyle.color}
                     onChange={(e) => setState("echart_features", index(), "itemStyle", "color", e.currentTarget.value)}/>
            </div>
          )
        }
        }
      </For>
    </div>);
}

function FeatureColoursInput() {
  return (
    <div class="mt-2">
      <label for="feature-colours">Feature colours</label>
      <HelpIcon
        helpMsg="Feature to colour map."/>
      <textarea id="feature-colours"
                class="w-full
                ml-1
                 mr-1
                 border
                 border-gray-300
                 dark:border-gray-700
                 rounded
                 min-h-[100px]
                 px-1
                 py-1
                 dark:bg-slate-900"
             value={
        state.echart_features.map(
          (feature) => {
            return `${feature.name}:${feature.itemStyle.color}`
        }).join("\n")
      }
             onChange={(e) => {
                const featureColours = e.currentTarget.value.split("\n").map((line) => {
                  const [feature, colour] = line.split(":");
                  return {feature, colour};
                });
                for (const featureColour of featureColours) {
                  const featureIndex = state.echart_features.findIndex((feature) => feature.name === featureColour.feature);
                  if (featureIndex >= 0) {
                    setState("echart_features", featureIndex, "itemStyle", "color", featureColour.colour);
                  }
                }
             }}
      />
    </div>
  )
}

export const ColourOptionsPanel: Component = () => {
  return <CollapsiblePanel title="Colour Options" defaultCollapsed={true}>
    <NuclColourSelect/>
    <FeatureColourPickers/>
    <FeatureColoursInput/>
  </CollapsiblePanel>
}