import {median, meanCoverage, medianCoverage, genomeCoverage} from "./coverageStat";
import {getCoverageChartOption} from "./echarts/getOption/getCoverageChartOption";
import {getGeneFeatureRenderer} from "./echarts/geneFeatures/getGeneFeatureRenderer";
import {getTooltips} from "./echarts/getOption/getTooltips";
// Import the echarts core module, which provides the necessary interfaces for using echarts.
import * as echarts from 'echarts/core';

import {LineChart, BarChart, CustomChart} from 'echarts/charts';
import {TooltipComponent, GridComponent, DataZoomComponent, DatasetComponent,} from 'echarts/components';
import {SVGRenderer, CanvasRenderer } from 'echarts/renderers';

echarts.use(
  [TooltipComponent, GridComponent,
      LineChart, BarChart,
      DataZoomComponent, CustomChart,
      DatasetComponent, SVGRenderer, CanvasRenderer]
);

export {getCoverageChartOption, median, meanCoverage, medianCoverage, genomeCoverage, getGeneFeatureRenderer, getTooltips, echarts};