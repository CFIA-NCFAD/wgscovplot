import {SampleDepths, VariantCall, WgsCovPlotDB} from "../db";
import {find, get, map} from "lodash";
import {state} from "../state";
import {ECColorArg} from "../db";

export function getVariantsSeries(db: WgsCovPlotDB) {
  const variantSeries = [];
  let i = 0;
  const depths = db.depths as SampleDepths;
  for (const sample of db.chartOptions.selectedSamples) {
    const sampleVariants = get(db.variants, sample, []) as VariantCall[];
    variantSeries.push({
      type: "bar",
      animation: false,
      xAxisIndex: i,
      yAxisIndex: i,
      data: map(sampleVariants, (x) => {
        const pos = typeof x.POS === "string" ? parseInt(x.POS) : x.POS;
        return [
          pos,
          depths[sample][pos - 1],
        ]
      }),
      barWidth: db.chartOptions.variantBarWidth,
      itemStyle: {
        color: function (arg: ECColorArg) {
          const pos = arg.data[0];
          const nt = find(sampleVariants, (x) => parseInt(x.POS) === pos)?.ALT;
          if (nt === undefined) {
            return "#333";
          }
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
        rotate: db.chartOptions.variantLabelsRotation,
        formatter: function (arg: ECColorArg) {
          const variant = find(sampleVariants, (x) => parseInt(x.POS) === arg.data[0]);
          if (variant === undefined) {
            return "";
          } else {
            const {REF, POS, ALT} = variant;
            return `${REF}${POS}${ALT}`;
          }
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