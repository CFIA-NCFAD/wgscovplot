import {whichSegment} from "./getSegmentsInfo";
import {find} from "lodash";
import {get, isNil} from "lodash";
import {SegmentCoords, VariantCall, WgsCovPlotDB} from "../../db";
import {ECColorArg} from "../../db";
import {state} from "../../state";

function segmentSeparatorLines(
  segments: string[] | undefined,
  segCoords: SegmentCoords | undefined
) {
  if (isNil(segments) || segments.length === 0) {
    return {};
  }
  if (isNil(segCoords)) {
    return {};
  }
  const data = [];
  for (let i = 0; i < segments.length; i++) {
    const segment = segments[i];
    const start = get(segCoords, [segment, "start"], Number.MAX_VALUE);
    const end = get(segCoords, [segment, "end"], Number.MAX_VALUE);
    if (start === Number.MAX_VALUE || end === Number.MAX_VALUE) {
      return {};
    }
    if (i === 0) {
      data.push({xAxis: end});
    } else if (i === segments.length - 1) {
      data.push({xAxis: start});
    } else {
      data.push({xAxis: start});
      data.push({xAxis: end});
    }
  }
  return {
    silent: true,
    symbol: ["none", "none"],
    label: {
      show: false,
    },
    lineStyle: {
      color: "#000",
      width: 1,
      type: "dashed",
      opacity: 0.2
    },
    data: data
  };
}

export function getSegmentVariantSeries(db: WgsCovPlotDB) {
  const variantSeries: object[] = [];
  let pos: number;
  const segCoords = db.segCoords;
  const segments = db.chartOptions.selectedSegments;
  if (isNil(segCoords) || isNil(segments) || segments.length === 0) {
    return variantSeries;
  }
  for (let i = 0; i < db.chartOptions.selectedSamples.length; i++) {
    const data = [];
    const sample = db.chartOptions.selectedSamples[i];
    for (const segment of segments) {
      const vars: undefined | VariantCall[]  = get(db, ["variants", sample, segment]);
      if (!isNil(vars)) {
        for (const varMap of vars.values()) {
          if (typeof varMap.POS === "string") {
            pos = parseInt(varMap.POS);
          } else {
            pos = varMap.POS
          }
          const depth = get(db, ["depths", sample, segment, pos - 1], 1E-7);
          const start = get(segCoords, [segment, 'start'], Number.MAX_VALUE);
          if (start === Number.MAX_VALUE) {
            continue;
          }
          data.push(
            [
              pos + start - 1, depth
            ]);
        }
      } else {
        data.push([]);
      }
    }
    const separatorLines = segmentSeparatorLines(db.chartOptions.selectedSegments, segCoords);
    variantSeries.push({
      type: "bar",
      xAxisIndex: i,
      yAxisIndex: i,
      data: data,
      animation: false,
      barWidth: db.chartOptions.variantBarWidth,
      itemStyle: {
        color: function (arg: ECColorArg) {
          const pos = arg.data[0];
          const segment = whichSegment(pos, db)
          const start = get(segCoords, [segment, 'start'], Number.MAX_VALUE);
          const seqPosition = pos - start + 1;
          const variants = get(db, ['variants', sample, segment], null);
          if (isNil(variants)) {
            return "#333";
          }
          const variant = find(variants, ["POS", seqPosition]);
          if (isNil(variant)) {
            return "#333";
          }
          const nt = variant.ALT_SEQ;
          return get(state.chartOptions.ntColor, nt, "#333");
        }
      },
      label: {
        show: db.chartOptions.showVariantLabels,
        position: "bottom",
        align: "left",
        verticalAlign: "middle",
        distance: 10,
        color: "inherit",
        rotate: db.chartOptions.variantLabelsRotation,
        formatter: function (arg: ECColorArg) {
          const pos = arg.data[0];
          const segment = whichSegment(pos, db);
          const start = get(segCoords, [segment, 'start'], Number.MAX_VALUE);
          const seqPosition = pos - start + 1;
          const variant = find(get(db.variants, [sample, segment]), {POS: seqPosition}, 0);
          if (!isNil(variant)) {
            const {REF_SEQ, POS, ALT_SEQ} = variant;
            return `${REF_SEQ}${POS}${ALT_SEQ}`;
          }
          return "";
        }
      },
      labelLayout: {
        hideOverlap: db.chartOptions.hideOverlappingVariantLabels
      },
      markLine: separatorLines,
      large: true,
      tooltip: {
        trigger: db.tooltipOptions.showTooltip ? "axis" : "none"
      }
    });
  }
  return variantSeries;
}