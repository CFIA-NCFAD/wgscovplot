import {median, meanCoverage, medianCoverage, genomeCoverage} from "./coverageStat";
import {getCoverageChartOption} from "./echarts/getOption/getCoverageChartOption";
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

export {getCoverageChartOption, median, meanCoverage, medianCoverage, genomeCoverage, echarts};