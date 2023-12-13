import {WgsCovPlotDB} from "../db";
import {isNil} from "lodash";

export const getDataZoom: (db: WgsCovPlotDB)
  => (any[] | [{ filterMode: string; xAxisIndex: number[]; type: string; startValue: number; endValue: number,  zoomLock: boolean }, { filterMode: string; xAxisIndex: number[]; showDataShadow: boolean; show: boolean; type: string; zoomLock: boolean }])
  = (db: WgsCovPlotDB) => {
  let xAxisIndex = [...Array(db.chartOptions.selectedSamples.length + 1).keys()];
  let start: number = 1
  let end: number = db.positions.length;
  if (!isNil(db.chart)) {
    if (!isNil(db.chart.getOption())) {
      let dataZoom: any = db.chart.getOption().dataZoom;
      start = Math.floor(dataZoom[0]["startValue"]);
      end = Math.floor(dataZoom[0]["endValue"]);
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