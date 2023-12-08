import {Component, Show} from "solid-js";
import {CollapsiblePanel} from "./CollapsiblePanel";
import {LowCovRegionsOptions} from "./LowCovRegionsOptions";
import {Slider} from "../components/Slider";
import {setState, state} from "../../state";
import {CovColour} from "./CovColour";
import {HelpIcon} from "../components/HelpIcon";
import {isNil} from "lodash";
import {SelectEChartsRenderer} from "./SelectEChartsRenderer";
import {ChartDarkModeToggle} from "./ChartDarkModeToggle";
import {NuclColourSelect} from "./NuclColourSelect";

export const DisplayOptionsPanel: Component = () => {
  return (
    <CollapsiblePanel title="Display Options" defaultCollapsed={true}>
      <LowCovRegionsOptions/>
      <Slider value={state.chartOptions.leftMargin}
              min={0} max={25} step={0.1}
              label="Left margin"
              helpInfo="Adjust the left margin of the plot."
              onChange={(e) => setState("chartOptions", "leftMargin", parseFloat(e.currentTarget.value))}/>
      <Slider value={state.chartOptions.rightMargin}
              min={0} max={25} step={0.1}
              label="Right margin"
              helpInfo="Adjust the right margin of the plot."
              onChange={(e) => setState("chartOptions", "rightMargin", parseFloat(e.currentTarget.value))}/>
      <Slider value={state.chartOptions.padTop}
              min={0} max={25} step={0.1}
              label="Top margin"
              helpInfo="Change the top margin of plot."
              onChange={(e) => setState("chartOptions", "padTop", parseFloat(e.currentTarget.value))}/>
      <Slider value={state.chartOptions.heightOffset}
              min={0} max={20} step={0.1}
              label="Vertical spacing"
              helpInfo="Adjust the vertical space between subplots if multiple samples are being shown."
              onChange={(e) => setState("chartOptions", "heightOffset", parseFloat(e.currentTarget.value))}/>
      <CovColour/>

      <div class="mt-2">
        <input type="checkbox" id="show-datazoom" class="hover:ring h-4 w-4 form-check"
               checked={state.chartOptions.showDataZoomSlider}
               onChange={(e) => setState("chartOptions", "showDataZoomSlider", e.currentTarget.checked)}/>
        <label class="ml-1" for="show-datazoom">Show data zoom slider?</label>
        <HelpIcon helpMsg="Show a slider to zoom in on a region of the genome."/>
      </div>
      <div class="mt-2">
        <label for="subplot-title-font-size">Subplot title font size</label>
        <HelpIcon helpMsg="Set the font size of the subplot titles."/>
        <input id="subplot-title-font-size" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
               type="number" min="0" step="0.5" value={state.chartOptions.subplotTitleFontSize}
               onChange={(e) => setState("chartOptions", "subplotTitleFontSize", parseFloat(e.currentTarget.value))}/>
      </div>
      <div class="mt-2">
        <label for="subplot-title-font-colour">Subplot title font colour</label>
        <HelpIcon helpMsg="Set the font colour of the subplot titles."/>
        <input id="subplot-title-font-colour"
               class="h-5 w-5 ml-2 hover:ring"
               type="color" value={state.chartOptions.subplotTitleColour}
               onChange={(e) => setState("chartOptions", "subplotTitleColour", e.currentTarget.value)}/>
      </div>
      <Show when={!isNil(state.amplicon_depths) && isNil(state.segments)}>
        <div class="mt-2">
          <input type="checkbox" id="show-amplicons" class="hover:ring h-4 w-4 form-check"
                 checked={state.show_amplicons}
                 onChange={(e) => setState("show_amplicons", e.currentTarget.checked)}/>
          <label class="ml-1" for="show-amplicons">Show Amplicon Depths</label>
        </div>
      </Show>
      <Show when={!isNil(state.primer_matches)}>
        <div class="mt-2">
          <input type="checkbox" id="show-primer-matches" class="hover:ring h-4 w-4 form-check"
                 checked={state.show_primer_matches}
                 onChange={(e) => setState("show_primer_matches", e.currentTarget.checked)}/>
          <label class="ml-1" for="show-primer-matches">Show Primer Matches</label>
        </div>
      </Show>
      <div class="mt-2">
        <input type="checkbox" id="show-features" class="hover:ring h-4 w-4 form-check"
               checked={state.chartOptions.showFeatures}
               onChange={(e) => setState("chartOptions", "showFeatures", e.currentTarget.checked)}/>
        <label class="ml-1" for="show-features">Show features</label>
      </div>
      <Show when={state.chartOptions.showFeatures}>
        <div class="mt-2">
          <input type="checkbox" id="show-gene-labels" class="hover:ring h-4 w-4 form-check"
                 checked={state.chartOptions.showGeneLabels}
                 onChange={(e) => setState("chartOptions", "showGeneLabels", e.currentTarget.checked)}/>
          <label class="ml-1" for="show-gene-labels">Show gene labels</label>
        </div>
        <Show when={state.chartOptions.showGeneLabels}>
          <div class="mt-2">
            <label for="gene-label-text-size">Gene label text size</label>
            <HelpIcon helpMsg="Set the size of the gene labels in the features plot."/>
            <input id="gene-label-text-size" class="w-1/6 ml-3 border border-gray-300 hover:ring rounded px-1"
                   type="number" min="0" step="0.5" value={state.chartOptions.geneLabelTextSize}
                   onChange={(e) => setState("chartOptions", "geneLabelTextSize", parseFloat(e.currentTarget.value))}/>
          </div>
        </Show>
        <div class="mt-2">
          <label for="feature-plot-height-scaling">Features subplot height scaling</label>
          <HelpIcon helpMsg="Adjust the height of features subplot."/>
          <input id="feature-plot-height-scaling" class="w-1/6 ml-3 border border-gray-300 rounded px-1"
                 type="number" min="0" step="1" value={state.chartOptions.featurePlotHeightScaling}
                 onChange={(e) => setState("chartOptions", "featurePlotHeightScaling", parseInt(e.currentTarget.value))}/>
          <span class="ml-1">%</span>
        </div>
      </Show>
      <SelectEChartsRenderer/>
      <ChartDarkModeToggle/>
      <NuclColourSelect/>
    </CollapsiblePanel>
  );
}