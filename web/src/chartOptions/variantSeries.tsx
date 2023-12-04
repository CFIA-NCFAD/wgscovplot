import {SampleDepths, VariantCall, WgsCovPlotDB} from "../db";
import {get, map} from "lodash";
import {state} from "../state";
import {ECColorArg} from "../db";

export function getVariantsSeries(db: WgsCovPlotDB) {
  console.log("Trigger getVariantsSeries")
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
        rotate: db.chartOptions.variantLabelsRotation,
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