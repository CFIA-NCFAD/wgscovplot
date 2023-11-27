import {whichSegment} from "./getSegmentsInfo";
import {isEmpty, find} from "lodash";
import {remove, map, get, isNil} from "lodash";
import {SegmentCoords, WgsCovPlotDB} from "../../db";
import {ECColorArg} from "../../chart";
import {state} from "../../state";

const segmentSeparatorLines = (segments: string[], segCoords: SegmentCoords): {} => {

  let data = [];
  for (let i = 0; i < segments.length; i++) {
    let segment = segments[i];
    if (i === 0) {
      data.push({xAxis: segCoords[segment].end});
    } else if (i === segments.length - 1) {
      data.push({xAxis: segCoords[segment].start});
    } else {
      data.push({xAxis: segCoords[segment].start});
      data.push({xAxis: segCoords[segment].end});
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
  let variantSeries = [];
  let pos;
  for (let i = 0; i < db.chartOptions.selectedSamples.length; i++) {
    let data = [];
    let sample = db.chartOptions.selectedSamples[i];
    for (let segment of db.chartOptions.selectedSegments) {
      // @ts-ignore
      let vars = db.variants[sample][segment];
      if (!isNil(vars)) {
        for (let [k, varMap] of vars.entries()) {
          pos = parseInt(varMap.POS);
          data.push(
            [
              // @ts-ignore
              pos + db.segCoords[segment].start - 1, db.depths[sample][segment][pos - 1]
            ]);
        }
      } else {
        data.push([]);
      }
    }
    variantSeries.push({
      type: "bar",
      xAxisIndex: i,
      yAxisIndex: i,
      data: data,
      barWidth: 2,
      itemStyle: {
        color: function (arg: ECColorArg) {
          let pos = arg.data[0];
          const segment = whichSegment(pos, db);
          const seqPosition = pos - db.segCoords[segment].start;
          const nt = db.segments_ref_seq[sample][segment][seqPosition];
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
        rotate: -30,
        formatter: function (arg: ECColorArg) {
          let output = "";
          let pos = arg.data[0];
          let segment = whichSegment(pos, db);
          let seqPosition = pos - db.segCoords[segment].start + 1;
          // @ts-ignore
          let variant = find(db.variants[sample][segment], {POS: seqPosition.toString()}, 0);
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
      markLine: segmentSeparatorLines(db.chartOptions.selectedSegments, db.segCoords),
      large: true,
      tooltip: {
        trigger: db.tooltipOptions.showTooltip ? "axis" : "none"
      }
    });
  }
  return variantSeries;
}