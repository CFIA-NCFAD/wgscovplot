import {median, meanCoverage, medianCoverage, genomeCoverage} from "./coverageStat";
import {getCoverageChartOption} from "./getOption/getCoverageChartOption";
import {getVariantHeatmapOption} from "./getOption/getVariantHeatmapOption";
import {getGeneFeatureRenderer} from "./geneFeatures/getGeneFeatureRenderer";
import {getRegionAmpliconDepthRenderer} from "./getOption/getRegionAmpliconDepthSeries";
import {getTooltips} from "./getOption/getTooltips";
import {getGrids} from "./getOption/getGrids";
import {getFluCoverageChartOption} from "./getOption/fluOption/getFluCoverageChartOption";
import {getMaxSegmentsLength} from "./getOption/fluOption/getFluSegmentsInfo";
import {getFluTooltips} from "./getOption/fluOption/getFluTooltips";
import "bootstrap/dist/css/bootstrap.css";
import "bootstrap/dist/js/bootstrap.js";
import $ from "jquery";
import JQuery from "jquery";
import {popper} from "@popperjs/core";
import select2 from "select2";
import "select2/dist/css/select2.css";
// Import the echarts core module, which provides the necessary interfaces for using echarts.
import * as echarts from "echarts/core";

import {LineChart, BarChart, CustomChart, HeatmapChart} from "echarts/charts";
import {
      TooltipComponent, GridComponent, DataZoomComponent,
      DatasetComponent, ToolboxComponent, VisualMapComponent, TitleComponent, MarkAreaComponent, MarkLineComponent
} from "echarts/components";
import {SVGRenderer, CanvasRenderer } from "echarts/renderers";

echarts.use(
  [TooltipComponent, GridComponent,
      LineChart, BarChart, ToolboxComponent, HeatmapChart,
      DataZoomComponent, CustomChart, VisualMapComponent, TitleComponent,
      DatasetComponent, SVGRenderer, CanvasRenderer, MarkAreaComponent, MarkLineComponent]
);

export {getCoverageChartOption, getFluCoverageChartOption, getVariantHeatmapOption,
      median, meanCoverage, medianCoverage, genomeCoverage,
      getGeneFeatureRenderer, getRegionAmpliconDepthRenderer,
      getTooltips, getGrids, getMaxSegmentsLength, getFluTooltips, echarts};