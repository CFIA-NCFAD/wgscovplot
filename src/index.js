import {median, meanCoverage, medianCoverage, genomeCoverage} from "./coverageStat";
import {getCoverageChartOption} from "./echarts/getOption/getCoverageChartOption";
import {getVariantHeatmapOption} from "./echarts/getOption/getVariantHeatmapOption";
import {getGeneFeatureRenderer} from "./echarts/geneFeatures/getGeneFeatureRenderer";
import {getRegionAmpliconDepthRenderer} from "./echarts/getOption/getRegionAmpliconDepthSeries";
import {getTooltips} from "./echarts/getOption/getTooltips";
import {getGrids} from "./echarts/getOption/getGrids";
import {getFluCoverageChartOption} from "./echarts/getOption/fluOptions/getFluCoverageChartOption";
import "bootstrap/dist/css/bootstrap.css";
import "bootstrap/dist/js/bootstrap.js";
import $ from "jquery";
import JQuery from "jquery";
import {popper} from "@popperjs/core";
import select2 from "select2";
import "select2/dist/css/select2.css";
// Import the echarts core module, which provides the necessary interfaces for using echarts.
import * as echarts from 'echarts/core';

import {LineChart, BarChart, CustomChart, HeatmapChart} from 'echarts/charts';
import {TooltipComponent, GridComponent, DataZoomComponent,
      DatasetComponent, ToolboxComponent, VisualMapComponent, TitleComponent} from 'echarts/components';
import {SVGRenderer, CanvasRenderer } from 'echarts/renderers';

echarts.use(
  [TooltipComponent, GridComponent,
      LineChart, BarChart, ToolboxComponent,HeatmapChart,
      DataZoomComponent, CustomChart,VisualMapComponent,TitleComponent,
      DatasetComponent, SVGRenderer, CanvasRenderer]
);

export {getCoverageChartOption, getFluCoverageChartOption, getVariantHeatmapOption, median, meanCoverage, medianCoverage, genomeCoverage,
      getGeneFeatureRenderer, getRegionAmpliconDepthRenderer, getTooltips, getGrids, echarts};