import {get, isNil, map} from "lodash";
import {getCoordsInterval} from "../util";
import {graphic} from "echarts/core";
import type {SampleDepths, SampleSegmentDepths, VariantCall, WgsCovPlotDB} from "../db";
import {SampleVariantCalls} from "../db";


/**
 * Renders the amplicon features shape
 * @param {boolean} show_amplicons - whether to plot amplicon feature or not (true for false)
 * @returns The shape of gene feature and/or amplicon feature
 */
export function getRegionAmpliconDepthRenderer(show_amplicons: boolean) {
  /**
   *
   * @param params - Echarts params
   * @param api - Echarts api
   */
  function renderRegionAmpliconDepth({coordSys}, api) {
    let [startX, startY] = api.coord([api.value(0), api.value(2)]);
    let [endX, endY] = api.coord([api.value(1), 1]);
    let rectShape = graphic.clipRectByRect(
      {
        x: startX,
        y: startY,
        width: endX - startX,
        height: endY - startY
      },
      coordSys,
    );
    return rectShape && {
      type: "rect",
      shape: rectShape,
      style: api.style(),
      invisible: !show_amplicons
    };
  }

  return renderRegionAmpliconDepth;
}

/**
 * Define options for amplicon depth coverage bars
 * @param {WgsCovPlotDB} db
 * @returns {Array<Object>}
 */
export const getRegionAmpliconDepthSeries = (db: WgsCovPlotDB): any[] => {
  let ampliconDepthSeries: any[] = [];
  if (isNil(db.amplicon_depths)) {
    return ampliconDepthSeries;
  }
  for (let [i, sample] of db.chartOptions.selectedSamples.entries()) {
    ampliconDepthSeries.push({
      type: "custom",
      xAxisIndex: i,
      yAxisIndex: i,
      renderItem: getRegionAmpliconDepthRenderer(show_amplicons),
      label: {
        show: false,
        position: "top",
        distance: 25,
        rotate: 60
      },
      labelLayout: {
        hideOverlap: false
      },
      encode: {
        x: [0, 1],
        y: 2,
      },
      tooltip: {
        trigger: "none"
      },
      silent: true,
      data: db.amplicon_depths[sample],
    });
  }
  return ampliconDepthSeries;
}
