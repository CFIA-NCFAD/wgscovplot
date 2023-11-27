import {getCustomXAxisLabel, getSegmentCoords} from "./chartOptions/segmented/getSegmentsInfo";
import {FEATURE_PLOT_PROPS, getCoordsInterval, hexToHSL, shapePoints, getTextWidth} from "./util";
import {genomeCoverage, meanCoverage, medianCoverage} from "./stats";
import {constant, find, get, isArray, isNil, map, times, sum, values, some} from "lodash";
import {getCoverageStatComparison, getVariantComparison} from "./chartOptions/tooltips";
import {ECFeature, ECFormatterFeature, SampleDepths, SegmentCoords, VariantCall, WgsCovPlotDB} from "./db";
import {toFloat16Array} from "./util";
import {setState, state} from "./state";
import {unwrap} from "solid-js/store";
import {graphic} from "echarts/core";
import {getFluGeneFeature, getMaxSegmentLength, whichSegment} from "./chartOptions/segmented/getSegmentsInfo";
import {getSegmentVariantComparison, getSegmentCoverageStatComparison} from "./chartOptions/tooltips";
import {getSegmentVariantSeries} from "./chartOptions/segmented/series";
import {getSegmentPrimerSeries, getPrimerInfo} from "./chartOptions/segmented/pcr";

type Grid = {
  show: boolean,
  height: string,
  top: string,
  left: string,
  right: string,
}


export interface ECColorArg {
  data: {
    0: number;
  }
}


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


export const tooltipFormatter = (db: WgsCovPlotDB) => {
  if (isNil(db.segments)) {
    let depths: SampleDepths = unwrap(db.depths) as SampleDepths;
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      let selectedSamples = db.chartOptions.selectedSamples;
      let output = "";
      let [{
        axisIndex,
        axisValue: position
      }] = params;
      if (axisIndex >= selectedSamples.length || db.variants === undefined) {
        return output;
      }
      let sample = selectedSamples[axisIndex];
      let sampleDepths = depths[sample];
      // change tooltip position if it's a click event and the tooltip is not already showing
      if (!db.tooltipOptions.show) {
        setState("tooltipOptions", "left", db.tooltipOptions.x + 10);
        setState("tooltipOptions", "top", db.tooltipOptions.y);
      }
      // set tooltip state
      setState("tooltipOptions", "show", true);
      setState("tooltipOptions", "sample", sample);
      setState("tooltipOptions", "position", position);
      setState("tooltipOptions", "depth", sampleDepths[position - 1]);

      let [dataZoom] = db.chart.getOption().dataZoom;
      let zoomStart = 1
      let zoomEnd = db.ref_seq.length;
      if (dataZoom !== undefined) {
        zoomStart = Math.floor(dataZoom.startValue);
        zoomEnd = Math.floor(dataZoom.endValue);
      }
      let positionRows: string[][] = [];
      let tables = [];
      const isVariantBar = params.find(x => x.componentSubType === "bar");
      if (isVariantBar) {
        if (db.chartOptions.crossSampleComparisonInTooltips) {
          positionRows = getVariantComparison(db, position);
        } else {
          let foundObj: any = find(db.variants[sample], {POS: position}, 0);
          if (!isNil(foundObj)) {
            for (const [key, value] of Object.entries(foundObj)) {
              positionRows.push([key, value as string]);
            }
          }
        }
        tables.push({headers: ["Variant Info", ""], rows: positionRows});
      } else {
        positionRows = [
          ["Sample", sample],
          ["Position", position],
          ["Sequence", db.ref_seq[position - 1]]
        ]
        tables.push({headers: ["Position Info", ""], rows: positionRows});
      }
      if (positionRows.length) {
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getCoverageStatComparison(db, zoomStart, zoomEnd, position);
          } else {
            let meanCov = meanCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            let medianCov = medianCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
            let genomeCov = genomeCoverage(sampleDepths, zoomStart, zoomEnd, db.chartOptions.low_coverage_threshold).toFixed(2);
            coverageStatRows = [
              [
                "Range",
                `${zoomStart.toLocaleString()} - ${zoomEnd.toLocaleString()}`,
              ],
              ["Mean Coverage", `${meanCov}X`],
              ["Median Coverage", `${medianCov}X`],
              [`Genome Coverage (>= ${db.chartOptions.low_coverage_threshold}X)`, `${genomeCov}%`],
            ];
          }
          tables.push({
            headers: ["Coverage View Stats", ""],
            rows: coverageStatRows,
          });
        }
      }
      setState("tooltipOptions", "tables", tables);
      return;
    };
  } else {
    return function (params: { axisIndex: number, axisValue: number, componentSubType: string }[]) {
      let [{
        axisIndex,
        axisValue: position
      }] = params;
      let positionRows: string[][] = [];
      let tables = [];
      let sample = db.chartOptions.selectedSamples[axisIndex];
      // find segment index which pos belongs to
      let segment = whichSegment(position, db);
      let sequence = db.segments_ref_seq[sample][segment];
      let segmentLength = sequence.length;
      // convert to pos in segment
      position = position - db.segCoords[segment].start + 1;
      let coverageDepth;
      let sampleSegDepths: number[] = db.depths[sample][segment];
      let refID = db.segments_ref_id[sample][segment];
      if (segmentLength === 0) {
        coverageDepth = `No result reported for segment ${segment}`;
      } else {
        if (position <= segmentLength) {
          // get coverage depth for pos
          coverageDepth = sampleSegDepths[position - 1].toLocaleString();
        } else {
          coverageDepth = `No sequence at this position. 
                        Reference sequence ${refID} 
                        is only ${sampleSegDepths.length} bp`;
        }
      }
      if (!db.tooltipOptions.show) {
        setState("tooltipOptions", "left", db.tooltipOptions.x + 10);
        setState("tooltipOptions", "top", db.tooltipOptions.y);
      }
      // set tooltip state
      setState("tooltipOptions", "show", true);
      setState("tooltipOptions", "sample", sample);
      setState("tooltipOptions", "position", position);
      setState("tooltipOptions", "depth", coverageDepth);
      const isVariantBar = params.find(element => element.componentSubType === "bar");
      if (isVariantBar) { //tooltips at Variant sites
        if (db.chartOptions.crossSampleComparisonInTooltips) {
          positionRows = getSegmentVariantComparison(db, sample, segment, position)
        } else {
          // @ts-ignore
          let foundObj = find(db.variants[sample][segment],
            {"POS": position.toString()});
          if (!isNil(foundObj)) {
            for (const [key, value] of Object.entries(foundObj)) {
              // @ts-ignore
              positionRows.push(...[[key, value]]);
            }
          }
          positionRows.push(["Coverage Depth", coverageDepth]);
        }
        let sortingCriteria = ["Sample", "Segment", "POS", "Segment Length", "Coverage Depth", "REF_ID", "REF_SEQ", "ALT_SEQ", "ALT_FREQ"]
        positionRows.sort((a, b) => sortingCriteria.indexOf(a[0]) - sortingCriteria.indexOf(b[0]))
        tables.push({headers: ["Variant Info", ""], rows: positionRows});
      } else { // tooltips for Non-Variant Sites
        if (position > segmentLength) { //Out of range when segment length < padding array
          positionRows = [
            ["Position", position],
            [coverageDepth, ""]
          ];
        } // Pos within segment length
        positionRows = [
          ["Sample", sample],
          ["Segment", segment],
          ["POS", position],
          ["Coverage Depth", coverageDepth],
          ["Segment Length", segmentLength],
          ["REF_ID", refID],
          ["REF_SEQ", sequence[position - 1]],
          ["ALT_SEQ", ""],
          ["ALT_FREQ", ""],
        ];
        tables.push({headers: ["Position Info", ""], rows: positionRows});
      }
      if (positionRows.length) { // write rows to table
        if (db.chartOptions.showCovStatsInTooltips) {
          let coverageStatRows = [];
          if (db.chartOptions.crossSampleComparisonInTooltips) {
            coverageStatRows = getSegmentCoverageStatComparison(db, segment, position);
          } else {
            let meanCov = meanCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            let medianCov = medianCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
            let genomeCov = genomeCoverage(sampleSegDepths, 1, segmentLength, db.chartOptions.low_coverage_threshold).toFixed(2);
            coverageStatRows = [
              ["Mean Coverage", `${meanCov}X`],
              ["Median Coverage", `${medianCov}X`],
              [`Genome Coverage (>= ${db.chartOptions.low_coverage_threshold}X)`, `${genomeCov}%`],
            ];
          }
          tables.push({headers: ["Coverage View Stats", ""], rows: coverageStatRows});
        }
        if (!isNil(db.primer_matches)) {
          let primerInfoRows = getPrimerInfo(sample, position, segment, db)
          if (primerInfoRows.length) {
            tables.push({headers: ["Primer Info", ""], rows: primerInfoRows});
          }
        }
      }
      setState("tooltipOptions", "tables", tables);
      return;
    }
  }
}

export const getTooltips = (db: WgsCovPlotDB) => {
  return [
    {
      trigger: "axis",
      enterable: true,
      triggerOn: db.chartOptions.tooltipTriggerOn,
      appendToBody: false,
      renderMode: "html",
      showContent: true,
      confine: true,
      position: "cursor",
      className: "hidden",
      axisPointer: {
        type: "line"
      },
      formatter: tooltipFormatter(db),
    },
  ];
}


export const getDatasets = (db: WgsCovPlotDB) => {
  //console.time("getDatasets");
  let datasets = [];
  if (!isNil(db.segments)) {     // segmented virus
    let maxSegmentLength = getMaxSegmentLength(state);
    setState("maxSegmentLength", maxSegmentLength);
    let positions = [...Array(sum(values(maxSegmentLength)) + 1).keys()];
    positions.shift()
    setState("positions", positions);
    let segCoords = getSegmentCoords(state);
    setState("segCoords", segCoords);
    let echart_features: ECFeature[] = getFluGeneFeature(state);
    setState("echart_features", echart_features)
    for (let sample of db.chartOptions.selectedSamples) {
      let depthArray: number [] = [];
      for (let segment of db.chartOptions.selectedSegments) {
        // @ts-ignore
        let ds = db.depths[sample][segment];
        let coords = db.segCoords[segment];
        if (ds.length < coords.maxLength) {
          // padding value 1E-5
          let padding = times(coords.maxLength - ds.length, constant(1E-5));
          ds = [...ds, ...padding]
        }
        depthArray = [...depthArray, ...ds];
      }
      datasets.push({
        dimensions: [
          {name: "depth", type: "float"},
          {name: "position", type: "int"},
        ],
        source: {
          position: db.positions,
          depth: depthArray,
        },
      });
    }
  } else {
    // non-segmented virus
    for (let sample of db.chartOptions.selectedSamples) {
      datasets.push({
        dimensions: [
          {name: "depth", type: "float"},
          {name: "position", type: "int"},
        ],
        source: {
          position: db.positions,
          depth: db.depths[sample],
        },
      });
    }
  }
  //console.log(datasets)
  //console.timeEnd("getDatasets");
  return datasets;
}

function getCoverageThresholdLine(db: WgsCovPlotDB) {
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

function getMarkArea(db: WgsCovPlotDB, sample: string) {
  let data = [];
  if (!(sample in db.depths)) {
    return {};
  }
  let threshold = db.chartOptions.low_coverage_threshold;
  // use unwrap Solid store util function to get underlying data for more rapid access to data
  let depths = unwrap(db.depths);
  if (!isNil(db.segments) && !isNil(db.segCoords) && db.segments.length > 0) {
    for (let segment of db.segments) {
      let sampleDepths = depths[sample];
      if (!(segment in sampleDepths) || isArray(sampleDepths)) {
        continue;
      }
      let sampleSegDepths: number[] = sampleDepths[segment] as number[];
      if (isNil(sampleSegDepths)) {
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
      show: false,
      position: "insideTop",
      fontSize: 10,
      rotate: 30,
      overflow: "truncate",
      ellipsis: "..."
    },
    data: data
  };
}

const getDepthSeries = (db: WgsCovPlotDB) => {
  //console.time("getDepthSeries")
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
  //console.timeEnd("getDepthSeries")
  return depthSeries;
}

export const getXAxes = (db: WgsCovPlotDB) => {
  let formatter: any = {};
  if (!isNil(db.segments) && db.segments.length > 0) {
    formatter.formatter = function (value: any) {
      return getCustomXAxisLabel(value, db.segments as string[], db.segCoords as SegmentCoords)
    }
  }
  let axes = [];
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
  let axes = [];
  for (let [i, sample] of db.chartOptions.selectedSamples.entries()) {
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
export const getGrids = (db: WgsCovPlotDB) => {
  let featureHeight: number = 0.0;
  if (db.chartOptions.showFeatures && db.show_amplicons && db.show_genes) {
    featureHeight = 15.0;
  } else if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    featureHeight = 6.0;
  } else {
    // no features subplot shown
    featureHeight = -5.0;
  }
  // TODO: find out what doubleStrand is
  featureHeight = (db.doubleStrand && featureHeight > 0) ? (featureHeight + 6.0) : featureHeight;
  featureHeight *= (db.chartOptions.featurePlotHeightScaling / 100);
  let padTop = db.chartOptions.padTop;
  let plotHeight = (db.chartOptions.showDataZoomSlider) ? 90 : 100;
  let subPlotHeight = plotHeight - featureHeight;
  let grids = [];
  let nSamples = db.chartOptions.selectedSamples.length;
  let heightOffset = db.chartOptions.heightOffset;
  let sampleHeight = (subPlotHeight - padTop) / nSamples - heightOffset;
  for (let idx = 0; idx < nSamples; idx++) {
    grids.push({
      show: true,
      height: sampleHeight + "%",
      top: ((sampleHeight + heightOffset) * idx + padTop) + "%",
      left: db.chartOptions.leftMargin + "%",
      right: db.chartOptions.rightMargin + "%",
    });
  }
  if (db.chartOptions.showFeatures && (db.show_amplicons || db.show_genes)) {
    grids.push({
      show: false,
      height: featureHeight + "%",
      top: ((sampleHeight + heightOffset) * nSamples + padTop) + "%",
      left: db.chartOptions.leftMargin + "%",
      right: db.chartOptions.rightMargin + "%",
    });
  }
  //console.log("grids", grids);
  return grids;
}

function getGeneFeatureRenderer(db: WgsCovPlotDB) {
  let yStart = 0;

  function renderGeneFeatures(params: { dataIndex: any; coordSys: any; }, api: {
    coord: (arg0: any[]) => [any, any] | [any];
    style: () => any;
  }) {
    const {
      dataIndex,
      coordSys,
    } = params;
    if (db.echart_features === undefined) {
      return null;
    }
    const feature: ECFeature = db.echart_features[dataIndex];

    const leftCoord = coordSys.x;
    const rightCoord = coordSys.width + coordSys.x;
    const [startX, startY] = api.coord([feature.value.start, dataIndex]);
    if (dataIndex === 0) {
      yStart = startY;
    }
    const [endX] = api.coord([feature.value.end, dataIndex]);
    const height = db.chartOptions.geneLabelTextSize + 3;
    const width = endX - startX;
    const y = yStart - height / 2 - feature.value.level;
    const points = shapePoints(startX, y, width, height, feature.value.strand, feature.value.type);
    if (points === null) {
      return null;
    }
    const textWidth = getTextWidth(feature.name, `normal ${db.chartOptions.geneLabelTextSize}px Arial`);
    const shape = graphic.clipPointsByRect(points, coordSys);
    let invisible = false;
    if (feature.value.type === "gene") {
      if (db.chartOptions.showGeneLabels) {
        // Element width is too small and hide label at the edges
        invisible = width < 10 || startX >= rightCoord || endX <= leftCoord;
      } else {
        invisible = true;
      }
      return {
        type: "polygon",
        shape: {
          points: shape,
        },
        style: api.style(),
        textContent: {
          type: "text",
          invisible: invisible,
          style: {
            text: feature.name,
            fill: width > textWidth ? (hexToHSL(feature.itemStyle.color).l > 50 ? "black" : "white") : "black",
            fontStyle: "normal",
            fontSize: db.chartOptions.geneLabelTextSize,
            fontWeight: "normal",
            stroke: feature.itemStyle.color,
            lineWidth: 1,
          },
        },
        textConfig: {
          position: width <= textWidth ? [0, -db.chartOptions.geneLabelTextSize / 2 - 1] : "inside",
          rotation: width <= textWidth ? 0.5 : 0,
          local: false,
        },
      };
    } else if (feature.value.type === "amplicon") {
      return {
        type: "polygon",
        shape: {
          points: shape,
        },
        style: api.style(),
        textContent: {},
        textConfig: {},
        invisible: !db.show_amplicons
      };
    } else if (feature.value.type === "segment") {
      invisible = width < 10 || startX >= rightCoord || endX <= leftCoord;
      return {
        type: "polygon",
        shape: {
          points: shape,
        },
        style: api.style(),
        textContent: {
          type: "text",
          invisible: invisible,
          style: {
            text: feature.name,
            fill: feature.itemStyle.color,
            fontStyle: "normal",
            fontSize: 12,
            fontWeight: "bolder",
          },
        },
        textConfig: {
          position: "top",
          distance: 18,
          rotation: db.chartOptions.geneLabelRotation,
          origin: "center",
          local: true,
        },
      };
    } else {
      return null;
    }
  }

  return renderGeneFeatures;
}

const getGeneFeatureSeries = (db: WgsCovPlotDB) => {
  let index = db.chartOptions.selectedSamples.length;
  return {
    type: "custom",
    animation: false,
    xAxisIndex: index,
    yAxisIndex: index,
    renderItem: getGeneFeatureRenderer(db),
    labelLayout: {
      hideOverlap: false,
    },
    data: db.echart_features,
    tooltip: {
      trigger: "item",
      show: db.tooltipOptions.showTooltip,
      enterable: true,
      appendToBody: true,
      triggerOn: "mousemove",
      renderMode: "html",
      borderRadius: 6,
      borderWidth: 2,
      showContent: "true",
      position: "top",
      textStyle: {
        fontSize: 14,
        fontWeight: "normal",
      },
      formatter: function (feature: ECFormatterFeature) {
        if (isNil(db.segments)) {
          return <div class="w-full h-full">
            <p class="text-sm font-bold">{/*@once*/ feature.name}</p>
            <p class="text-sm">
              {/*@once*/ feature.value.start.toLocaleString()} - {/*@once*/ feature.value.end.toLocaleString()}
            </p>
            <p class="text-sm">{/*@once*/ (feature.value.end - feature.value.start + 1).toLocaleString()} bp</p>
            <button class="hover:ring py-1 px-2 rounded mt-2"
                    style={{
                      "background-color": /*@once*/ feature.color !== undefined ? feature.color : "#2a2a2a",
                      "color": /*@once*/ feature.color !== undefined ? (hexToHSL(feature.color).l > 50 ? "black" : "white") : "#fff",
                    }}
                    onClick={(e) => {
                      const seq = db.ref_seq.slice(feature.value.start - 1, feature.value.end);
                      // TODO: add reference ID/name to header; need to make sure that it is exported from the backend
                      const header = `>REFID|${feature.value.start}-${feature.value.end} ${feature.name}`;
                      navigator.clipboard.writeText(`${header}\n${seq}\n`);
                      // Indicate that the text has been copied
                      // Using a Solid signal will produce the following warning:
                      // "computations created outside a `createRoot` or `render` will never be disposed"
                      // so the button text is updated directly
                      e.currentTarget.innerHTML = "Copied seq!";
                    }}>
              Copy seq
            </button>
          </div>;
        }
      },
    },
  };
}

export function getVariantsSeries(db: WgsCovPlotDB) {
  let variantSeries = [];
  let i = 0;
  const depths = db.depths as SampleDepths;
  for (let sample of db.chartOptions.selectedSamples) {
    let sampleVariants = get(db.variants, sample, []) as VariantCall[];
    variantSeries.push({
      type: "bar",
      animation: false,
      xAxisIndex: i,
      yAxisIndex: i,
      data: map(sampleVariants, (x) => {
        let pos = parseInt(x.POS);
        return [
          pos,
          depths[sample][pos - 1],
        ]
      }),
      barWidth: 2,
      itemStyle: {
        color: function (arg: ECColorArg) {
          let pos = arg.data[0];
          let nt = db.ref_seq[pos - 1];
          return get(state.chartOptions.ntColor, nt, "#333");
        },
      },
      label: {
        show: db.chartOptions.showVariantLabels,
        position: "bottom",
        align: "left",
        verticalAlign: "middle",
        distance: 10,
        color: "inherit",
        rotate: -30,
        formatter: function (arg: ECColorArg) {
          let output = "";
          Object.values(sampleVariants).forEach(({POS, REF, ALT}) => {
            let pos = `${arg.data[0]}`;
            if (parseInt(POS) === parseInt(pos)) {
              output += `${REF}${POS}${ALT}`;
            }
          });
          return output;
        }
      },
      labelLayout: {
        hideOverlap: db.chartOptions.hideOverlappingVariantLabels
      },
      tooltip: {
        trigger: db.tooltipOptions.showTooltip ? "axis" : "none"
      }
    });
    i++;
  }
  return variantSeries;
}

export function getRegionAmpliconDepthRenderer(db: WgsCovPlotDB) {
  // @ts-ignore
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
      invisible: !db.show_amplicons
    };
  }

  return renderRegionAmpliconDepth;
}

export function getRegionAmpliconDepthSeries(db: WgsCovPlotDB) {
  let ampliconDepthSeries: any[] = [];
  if (isNil(db.amplicon_depths))
    return ampliconDepthSeries;
  for (let [i, sample] of db.chartOptions.selectedSamples.entries()) {
    ampliconDepthSeries.push({
      type: "custom",
      xAxisIndex: i,
      yAxisIndex: i,
      renderItem: getRegionAmpliconDepthRenderer(db),
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

export const getSeries = (db: WgsCovPlotDB) => {
  let series: any[] = [];
  series.push(...getDepthSeries(db));
  if (!isNil(db.amplicon_depths)) {
    series.push(...getRegionAmpliconDepthSeries(db));
  }
  if (!isNil(db.variants) && db.chartOptions.showVariants) {
    if (isNil(db.segments)) {
      series.push(...getVariantsSeries(db));
    } else {
      series.push(...getSegmentVariantSeries(db));
    }
  }
  if (!isNil(db.primer_matches)) {
    series.push(...getSegmentPrimerSeries(db));
  }
  let isDoubleStrand: boolean = some(db.echart_features, {value: {strand: -1}})
  setState("doubleStrand", isDoubleStrand)
  if (db.chartOptions.showFeatures && (db.show_genes || db.show_amplicons)) {
    let geneFeatureSeries: any = getGeneFeatureSeries(db);
    series.push(geneFeatureSeries);
  }
  return series;
}