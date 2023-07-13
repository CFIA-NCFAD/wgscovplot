import {Component} from "solid-js";
import {HelpIcon} from "../components/HelpIcon";
import {setState, state} from "../../state";

export const CovColour: Component = () => {
  return <div class="w-full mt-2 inline-block">
    <label for="cov-colour" class="">Coverage colour</label>
    <HelpIcon helpMsg="Set the main colour of coverage plots."/>
    <input type="color" value={state.chartOptions.covColour}
           class="hover:ring h-5 w-5" id="cov-colour"
           onChange={(e) => setState("chartOptions", "covColour", e.currentTarget.value)}/>
  </div>
}