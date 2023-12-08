import {WgsCovPlotDB} from "../db";

export const getDataZoom: (db: WgsCovPlotDB)
  => (any[] | [{ filterMode: string; xAxisIndex: number[]; type: string; zoomLock: boolean }, { filterMode: string; xAxisIndex: number[]; showDataShadow: boolean; show: boolean; type: string; zoomLock: boolean }])
  = (db: WgsCovPlotDB) => {
  let xAxisIndex = [...Array(db.chartOptions.selectedSamples.length + 1).keys()];
  return [
    {
      type: "inside",
      filterMode: "none",
      xAxisIndex: xAxisIndex,
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