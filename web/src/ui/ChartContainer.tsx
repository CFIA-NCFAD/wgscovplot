import {setState, state} from "../state";
import type {Component} from 'solid-js';
import {createEffect, createMemo, onMount, untrack} from "solid-js";
import * as echarts from "echarts/core";

import {
  DatasetComponent,
  DataZoomComponent,
  GridComponent,
  MarkAreaComponent,
  MarkLineComponent,
  TitleComponent,
  ToolboxComponent,
  TooltipComponent,
  VisualMapComponent,
} from "echarts/components";
import {LabelLayout} from "echarts/features";
import {BarChart, CustomChart, HeatmapChart, LineChart} from "echarts/charts";
import {CanvasRenderer, SVGRenderer} from "echarts/renderers";
import type {ECharts, EChartsType} from "echarts";
import {isNil, map, max, pick, sum, values} from "lodash";
import {MosdepthInfo} from "../db";
import {unwrap} from "solid-js/store";
import {getSeries} from "../chartOptions/chartSeries";
import {getXAxes, getYAxes} from "../chartOptions/axes";
import {getGrids} from "../chartOptions/chartGrid";
import {getDataZoom} from "../chartOptions/dataZoom";
import {getDatasets} from "../chartOptions/dataset";
import {getTooltips} from "../chartOptions/tooltips";
import {DragDropProvider, DragDropSensors, DragEventHandler,} from "@thisbeyond/solid-dnd";
import {ChartTooltip} from "./ChartTooltip";


echarts.use(
  [TooltipComponent, GridComponent,
    LineChart, BarChart, ToolboxComponent, HeatmapChart,
    DataZoomComponent, CustomChart, VisualMapComponent, TitleComponent,
    DatasetComponent, SVGRenderer, CanvasRenderer, MarkAreaComponent,
    MarkLineComponent, LabelLayout]
);

// add global event listeners for mouse position
document.addEventListener('mousemove', onMouseUpdate, false);
document.addEventListener('mouseenter', onMouseUpdate, false);

function onMouseUpdate(e: MouseEvent) {
  setState("tooltipOptions", "x", e.pageX);
  setState("tooltipOptions", "y", e.pageY);
}

// set state.chartOptions.yMax (the max y-axis value) based on the max depth of the selected samples (and segments)
createEffect(() => {
  if (isNil(state.mosdepth_info) || isNil(state.chartOptions.selectedSamples)) {
    return;
  }
  let sampleMaxDepth = map(state.chartOptions.selectedSamples, (sample) => {
    let info: { [key: string]: MosdepthInfo } | MosdepthInfo = state.mosdepth_info[sample];
    if (isNil(info)) {
      return 0;
    }
    // for segmented viruses, we need to get the max depth of all selected segments
    if (!isNil(state.chartOptions.selectedSegments) && !isNil(state.segCoords)) {
      let infos: MosdepthInfo[] = values(pick(info as { [key: string]: MosdepthInfo }, state.chartOptions.selectedSegments));
      return max(map(infos, (i) => i.max_depth));
    }
    // for non-segmented viruses, we can just get the max depth of the sample
    return info.max_depth
  });
  let ymax = max(values(sampleMaxDepth)) as number;
  setState("chartOptions", "yMax", ymax);
});

createEffect(() => {
  if (isNil(state.segments)) return;
  let selectedSegments = state.chartOptions.selectedSegments;
  let selectedSamples = state.chartOptions.selectedSamples;
  untrack(() => {
    if (isNil(selectedSegments) || selectedSegments.length === 0) {
      setState("chartOptions", "selectedSegments", state.segments);
    }
    if (isNil(selectedSamples) || selectedSamples.length === 0) {
      setState("chartOptions", "selectedSamples", state.samples.slice(0, 3));
    }
  });
});

const Chart: Component = () => {
  let chartDiv: undefined | HTMLDivElement = undefined;
  let chart: undefined | ECharts | EChartsType = undefined;


  const initChart = () => {
    if (chartDiv !== undefined) {
      // if the chart has already been initialized, dispose of it so that it can be re-initialized
      if (chart !== undefined) {
        // @ts-ignore
        echarts.dispose(chart);
      }
      // @ts-ignore
      chart = echarts.init(chartDiv, state.chartOptions.darkMode ? "dark" : "white", {renderer: state.chartOptions.renderer});
      if (chart === undefined) {
        throw new Error("chart is undefined")
      }
      // @ts-ignore
      chart.on("click", function ({componentIndex, componentSubType, value: {start, end}}) {
        // @ts-ignore
        if (componentIndex === chart.getOption().series.length - 1 && componentSubType === "custom") {
          setState("chartOptions", "startPos", start);
          setState("chartOptions", "endPos", end);
        }
      });

      chart.on("dblclick", function ({componentIndex, componentSubType}) {
        // @ts-ignore
        if (componentIndex === chart.getOption().series.length - 1 && componentSubType === "custom") {
          const start = 1;
          const end = state.positions.length;
          setState("chartOptions", "startPos", start);
          setState("chartOptions", "endPos", end);
        }
      });
      setState("chart", chart);
    }
  }

  onMount(() => {
    initChart();
    // use built-in ResizeObserver to resize chart on div resize (https://developer.mozilla.org/en-US/docs/Web/API/ResizeObserver)
    const resizeObserver = new ResizeObserver(() => {
      if (chart !== undefined) {
        chart.resize();
      }
    });
    if (chartDiv !== undefined) {
      resizeObserver.observe(chartDiv);
    }
  });

  // re-init chart if renderer changes or if dark mode toggled
  createEffect(() => {
    state.chartOptions.renderer;
    state.chartOptions.darkMode;
    // don't track any other state changes
    untrack(() => {
      initChart();
      if (chart !== undefined) {
        chart.setOption(unwrap(echartsOptions()), {
          replaceMerge: ['dataset', 'xAxis', 'yAxis', 'series', 'grid', 'dataZoom']
        });
      }
    });
  });


  let datasets = createMemo(() => {
    return getDatasets(state);
  });
  let xAxes = createMemo(() => {
    return getXAxes(state);
  });
  let yAxes = createMemo(() => {
    return getYAxes(state);
  });
  let grids = createMemo(() => {
    return getGrids(state);
  });
  let series = createMemo(() => {
    return getSeries(state);
  });

  let echartsOptions = createMemo(() => {
    //console.time("build opts")
    let opts = {
      dataset: datasets(),
      xAxis: xAxes(),
      yAxis: yAxes(),
      series: series(),
      grid: grids(),
      tooltip: getTooltips(state),
      toolbox: {
        show: "true",
        feature: {
          saveAsImage: {
            name: "wgscovplot",
          },
        },
      },
      dataZoom: getDataZoom(state),
    };
    //console.timeEnd("build opts")
    return opts;
  });

  createEffect(() => {
    if (chart !== undefined) {
      //console.time("chart setOption")
      let dataZoom : any = chart.getOption().dataZoom;
      let start = Math.floor(dataZoom[0]["startValue"]);
      let end = Math.floor(dataZoom[0]["endValue"]);
      // tell ECharts to merge replace rather than normal merge
      // unwrap opts to pass a plain object instead of a Solid Proxy object to ECharts
      // this reduces deep cloning by ECharts and speeds up chart updates
      chart.setOption(unwrap(echartsOptions()), {
        replaceMerge: ['dataset', 'xAxis', 'yAxis', 'series', 'grid', 'dataZoom']
      });

      state.chart.dispatchAction({
        type: "dataZoom",
        startValue: start,
        endValue: end,
      });

      //console.log("Chart Options", chart.getOption())
      //console.timeEnd("chart setOption")
     // console.info("CHART UPDATE", Date.now());
    }
  })

  createEffect(() => {
    if (chart !== undefined) {
      state.chartOptions.sidebarCollapsed;
      chart.resize();
    }
  })

  // update dataZoom when startPos or endPos changes
  // dataZoom needs to be set in the chart options, i.e. state.chartOptions.showDataZoomSlider must be true
  createEffect(() => {
    if (state.chart !== undefined && state.chartOptions.startPos !== 0 && state.chartOptions.endPos !== 0) {
      let start = state.chartOptions.startPos;
      let end = state.chartOptions.endPos;
      state.chart.dispatchAction({
        type: "dataZoom",
        startValue: start,
        endValue: end,
      });
    }
  });

  return <div ref={chartDiv} class="flex-grow flex-shrink min-w-fit"></div>
}

export const ChartContainer: Component = () => {
  let transform = {x: 0, y: 0};
  const onDragMove: DragEventHandler = (e) => {
    transform = {...e.draggable.transform}
  };

  const onDragEnd: DragEventHandler = ({draggable}) => {
    const node = draggable.node;
    setState("tooltipOptions", "top", node.offsetTop + transform.y);
    setState("tooltipOptions", "left", node.offsetLeft + transform.x);
  };

  return <DragDropProvider onDragMove={onDragMove} onDragEnd={onDragEnd}>
    <DragDropSensors/>
    <main class={"flex " + (state.chartOptions.sidebarCollapsed ? "p-4 w-screen" : "sm:w-4/5 p-4")}>
      <Chart/>
      <ChartTooltip/>
    </main>
  </DragDropProvider>
}