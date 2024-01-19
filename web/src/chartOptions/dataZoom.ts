import {WgsCovPlotDB} from "../db";
import {get, isNil} from "lodash";

export function getDataZoom(db: WgsCovPlotDB) {
  const xAxisIndex = [...Array(db.chartOptions.selectedSamples.length + 1).keys()];
  let start: number = 1
  let end: number = db.positions.length;
  if (!isNil(db.chart)) {
    const opts = db.chart.getOption();
    const dataZoomStart = get(opts, ["dataZoom", 0, "startValue"], start);
    const dataZoomEnd = get(opts, ["dataZoom", 0, "endValue"], end);
    if (dataZoomEnd - dataZoomStart >= 5) {
      start = dataZoomStart;
      end = dataZoomEnd;
    }
  }
  return [
    {
      type: "inside",
      filterMode: "none",
      xAxisIndex: xAxisIndex,
      startValue: start,
      endValue: end,
      zoomLock: false,
    },
    {
      show: db.chartOptions.showDataZoomSlider,
      filterMode: "none",
      xAxisIndex: xAxisIndex,
      type: "slider",
      zoomLock: false,
      showDataShadow: false,
      showDetail: true,
    },
  ];
}