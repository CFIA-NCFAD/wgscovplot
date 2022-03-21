import {median, meanCoverage, medianCoverage, genomeCoverage} from "./coverageStat";
import {getCoverageChartOption} from "./echarts/getOption/getCoverageChartOption";
import {getVariantHeatmapOption} from "./echarts/getOption/getVariantHeatmapOption";
import {getGeneFeatureRenderer} from "./echarts/geneFeatures/getGeneFeatureRenderer";
import {getTooltips} from "./echarts/getOption/getTooltips";
import {getGrids} from "./echarts/getOption/getGrids";
// Import the echarts core module, which provides the necessary interfaces for using echarts.
import * as echarts from 'echarts/core';


import {LineChart, BarChart, CustomChart, HeatmapChart} from 'echarts/charts';
import {TooltipComponent, GridComponent, DataZoomComponent,
      DatasetComponent, ToolboxComponent, VisualMapComponent} from 'echarts/components';
import {SVGRenderer, CanvasRenderer } from 'echarts/renderers';

echarts.use(
  [TooltipComponent, GridComponent,
      LineChart, BarChart, ToolboxComponent,HeatmapChart,
      DataZoomComponent, CustomChart,VisualMapComponent,
      DatasetComponent, SVGRenderer, CanvasRenderer]
);

export {getCoverageChartOption, getVariantHeatmapOption, median, meanCoverage, medianCoverage, genomeCoverage,
      getGeneFeatureRenderer, getTooltips, getGrids, echarts};