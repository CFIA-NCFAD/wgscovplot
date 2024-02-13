import {Component} from "solid-js";
import {HelpIcon} from "../components/HelpIcon";
import {setState, state} from "../../state";

export const LowCovColour: Component = () => {
  return <div class="w-full mt-2 inline-block">
    <label for="low-cov-colour" class="">Colour</label>
    <HelpIcon helpMsg="Set the colour of low coverage regions in the coverage plots."/>
    <input type="color" value={state.chartOptions.lowCovColour}
           class="hover:ring h-5 w-5" id="low-cov-colour"
           onChange={(e) => setState("chartOptions", "lowCovColour", e.currentTarget.value)}/>
  </div>
}