import {WgsCovPlotDB} from "../db";
import {unwrap} from "solid-js/store";
import {isArray, isEmpty, isNil} from "lodash";
import {getCoordsInterval} from "../util";


export function getCoverageThresholdLine(db: WgsCovPlotDB) {
  return {
    silent: true,
    symbol: ["none", "none"],
    label: {
      show: true,
      formatter: "{c}X",
    },
    lineStyle: {
      color: db.chartOptions.lowCovThresholdLineColour,
      width: db.chartOptions.lowCovThresholdLineWidth,
      type: "dotted",
      opacity: 1
    },
    data: [
      {
        name: "Low Coverage Threshold",
        yAxis: db.chartOptions.low_coverage_threshold
      }
    ]
  };
}

export function getMarkArea(db: WgsCovPlotDB, sample: string) {
  let data = [];
  if (!(sample in db.depths)) {
    return {};
  }
  let threshold = db.chartOptions.low_coverage_threshold;
  // use unwrap Solid store util function to get underlying data for more rapid access to data
  let depths = unwrap(db.depths);
  if (!isNil(db.segments) && !isEmpty(db.segCoords) && db.segments.length > 0) {
    let segments = Object.keys(db.segCoords)
    for (let segment of segments) {
      let sampleDepths = depths[sample];
      if (!(segment in sampleDepths) || isArray(sampleDepths)) {
        continue;
      }
      let sampleSegDepths: number[] = sampleDepths[segment] as number[];
      if (isEmpty(sampleSegDepths)) {
        continue;
      }
      for (let [start, end] of getCoordsInterval(sampleSegDepths, threshold)) {
        data.push([
          {
            name: `${start}-${end} (<${db.chartOptions.low_coverage_threshold}X)`,
            xAxis: start + db.segCoords[segment].start - 1,
          },
          {
            xAxis: end + db.segCoords[segment].start - 1,
          }
        ]);
      }
    }
  } else {
    let sampleDepths: number[] = depths[sample] as number[];
    let intervals = getCoordsInterval(sampleDepths, threshold);
    for (let [start, end] of intervals) {
      data.push([
        {
          name: `${start}-${end} (<${db.chartOptions.low_coverage_threshold}X)`,
          xAxis: start,
        },
        {
          xAxis: end,
        }
      ]);
    }
  }
  return {
    itemStyle: {
      color: db.chartOptions.lowCovColour,
      opacity: db.chartOptions.lowCoverageOpacity,
    },
    label: {
      show: db.chartOptions.showLowCoverageCoords,
      position: "insideTop",
      fontSize: 10,
      rotate: db.chartOptions.coordsLabelsRotation,
      overflow: "truncate",
      ellipsis: "..."
    },
    data: data
  };
}

export const getDepthSeries = (db: WgsCovPlotDB) => {
  let depthSeries = [];
  for (let i = 0; i < db.chartOptions.selectedSamples.length; i++) {
    let sample = db.chartOptions.selectedSamples[i];
    let series = {
      type: "line",
      xAxisIndex: i,
      yAxisIndex: i,
      areaStyle: {
        color: db.chartOptions.covColour,
      },
      encode: {
        x: "position",
        y: "depth",
      },
      symbol: "none",
      datasetIndex: i,
      lineStyle: {
        color: db.chartOptions.covColour,
        opacity: 0,
      },
      tooltip: {
        trigger: db.tooltipOptions.variantSitesOnly ? "none" : "axis",
      },
      silent: true,
      large: true,
    };
    if (db.chartOptions.showLowCovRegions) {
      series = {
        ...{markArea: getMarkArea(db, sample)},
        ...series,
        ...{markLine: getCoverageThresholdLine(db),}
      }
    }
    depthSeries.push(series);
  }
  return depthSeries;
}