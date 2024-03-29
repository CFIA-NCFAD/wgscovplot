import {SegmentCoords, WgsCovPlotDB} from "../db";
import {isNil} from "lodash";
import {getCustomXAxisLabel} from "./segmented/getSegmentsInfo";
import {FEATURE_PLOT_PROPS} from "../util";


export const getXAxes = (db: WgsCovPlotDB) => {
  // eslint-disable-next-line
  const formatter: any = {};
  if (!isNil(db.segCoords)) {
    const segments = Object.keys(db.segCoords)
    if (segments.length > 0) {
      // eslint-disable-next-line
      formatter.formatter = function (value: any) {
        return getCustomXAxisLabel(value, segments, db.segCoords as SegmentCoords)
      }
    }
  }
  const axes = [];
  for (let i = 0; i < db.chartOptions.selectedSamples.length; i++) {
    axes.push({
      type: "value",
      gridIndex: i,
      min: 1,
      max: db.positions.length,
      minorTick: {show: true},
      axisLabel: {
        show: db.chartOptions.showXAxisLabel,
        interval: "auto",
        ...formatter,
      }
    });
  }
  if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    axes.push({
      type: "value",
      gridIndex: db.chartOptions.selectedSamples.length,
      min: 1,
      max: db.positions.length,
      axisLabel: {
        interval: "auto",
        ...formatter,
      },
    });
  }
  return axes;
}

export const getYAxes = (db: WgsCovPlotDB) => {
  const axes = [];
  for (const [i, sample] of db.chartOptions.selectedSamples.entries()) {
    axes.push({
      type: db.chartOptions.scaleType,
      gridIndex: i,
      name: sample,
      nameTextStyle: {
        fontStyle: "normal",
        fontWeight: "bolder",
        fontSize: db.chartOptions.subplotTitleFontSize,
        color: db.chartOptions.subplotTitleColour,
      },
      nameLocation: "end",
      nameRotate: 0.01,
      min: db.chartOptions.scaleType === "log" ? 1 : 0,
      max: db.chartOptions.yMax === 0 ? 10000 : db.chartOptions.yMax,
      minorSplitLine: {
        show: true,
      },
    });
  }
  if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    axes.push({
      max: FEATURE_PLOT_PROPS.max_grid_height,
      gridIndex: db.chartOptions.selectedSamples.length,
      show: false,
    });
  }
  return axes;
}