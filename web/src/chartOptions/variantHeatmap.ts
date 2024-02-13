import {toTableHtml} from "../util";
import {filter, find, isEmpty, isNil, map, orderBy, union, uniq} from "lodash";

import {WgsCovPlotDB} from "../db";

// eslint-disable-next-line
function getTooltipHeatmap(dataSamples: string[], mutations: any[], db: WgsCovPlotDB) {
  return {
    enterable: true,
    appendToBody: true,
    renderMode: "html",
    showContent: true,
    confine: true,
    formatter: function (arg: { value: never[] | undefined }) {
      let output = "";
      let mutation = "";
      let sample = "";
      // values may be undefined when mouse move and zoom/in out heat map (causing error log)
      const value = arg.value;
      if (!isNil(value)) {
        mutation = mutations[value[0]];
        sample = dataSamples[value[1]];
      }
      let rows: string[][] = [];
      const sampleVariants = db.variants[sample];
      if (!isNil(sampleVariants)) {
        const variantOfInterest: object | undefined = find(Object.values(sampleVariants), {mutation});
        if (!isNil(variantOfInterest)) {
          rows = [...Object.entries(variantOfInterest)];
        }
        if (!isEmpty(rows)) {
          output += toTableHtml({
            headers: ["", ""],
            rows,
          });
        } else {
          output += `Mutation ${mutation.bold()} was not found in sample ${sample.bold()}`;
        }
      } else {
        output += `No variant calling results for ${sample.bold()}`;
      }
      return output;
    }
  }
}

function getMutationMatrix(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  let sampleVariants: any[] = [];
  db.chartOptions.selectedSamples.forEach(sample => {
    const variantCalls = db.variants[sample]
    sampleVariants = union(sampleVariants, filter(variantCalls, "mutation"));
  });
  sampleVariants = orderBy(sampleVariants, "POS", "asc");
  const mutations = uniq(map(sampleVariants, "mutation"));
  const altFreqMatrix = [];
  for (let i = 0; i < db.chartOptions.selectedSamples.length; i++) {
    const sample = db.chartOptions.selectedSamples[i];
    for (let j = 0; j < mutations.length; j++) {
      const mutation = mutations[j];
      const foundObj = find(sampleVariants, {sample, mutation});
      altFreqMatrix.push([
        j,
        db.chartOptions.selectedSamples.length - 1 - i,
        !isNil(foundObj) ? foundObj.ALT_FREQ : 0
      ]);
    }
  }
  return [mutations, altFreqMatrix];
}

/**
 * Define all options for variant heatmap
 * @param {WgsCovPlotDB} db
 * @returns {Object} - The options for variant heatmap
 *
 * samples = ["sample1", "sample2","sample3", ..... "sampleN] // Dim Y
 * mutations = ["C19T", "T847C(orf1ab:D194D)","T3442C(orf1ab:N1059N)", ... , "mutation_n_name"] // Dim X
 * Please see https://echarts.apache.org/en/option.html#series-heatmap.data
 * 2D matrix in which mutations represent columns, samples represent rows
 * * variantMatrix = [
 *      //[dimX, dimY, heatmap value]
 *       [0, 1, 0.45],
 *       [0, 2, 0.45],
 *       [1, 0, 1]
 *   ]
 * * variants = {SRR17230678: Array(50), SRR17230672: Array(54), SRR17230675: Array(56), SRR17230680: Array(43), SRR17230671: Array(49)}
 SRR17230671: Array(49)
 0: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G210T', POS: 210, REF: 'G', …}
 1: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C241T', POS: 241, REF: 'C', …}
 2: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C3037T(orf1ab:F924F)', POS: 3037, REF: 'C', …}
 */
export function getVariantHeatmapOption(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  const [mutations, matrix]: any[][] = getMutationMatrix(db);
  const dataSamples = [...db.chartOptions.selectedSamples].reverse();
  return {
    xAxis: {
      type: "category",
      data: mutations,
      splitArea: {
        show: true
      },
      position: "bottom",
      axisLabel: {
        rotate: -45,
      },
    },
    yAxis: {
      type: "category",
      data: dataSamples,
      splitArea: {
        show: true
      },
    },
    dataZoom: [
      {
        type: "inside"
      },
      {
        type: "slider",
        show: db.chartOptions.showDataZoomSlider,
      }
    ],
    visualMap: {
      min: 0.00,
      max: 1.00,
      precision: 2,
      calculable: true,
      orient: "vertical",
      left: "right",
      top: "5%",
      inRange: {
        color: [
          "#ffffe5",
          "#f7fcb9",
          "#d9f0a3",
          "#addd8e",
          "#78c679",
          "#41ab5d",
          "#238443",
          "#005a32",
        ]
      }
    },
    series: [
      {
        type: "heatmap",
        data: matrix,
        label: {
          show: false
        },
        emphasis: {
          itemStyle: {
            shadowBlur: 10,
            shadowColor: "rgba(0, 0, 0, 0.5)"
          }
        }
      }
    ],
    tooltip: getTooltipHeatmap(dataSamples, mutations, db),
    toolbox: {
      show: "true",
      feature: {
        dataView: {
          readOnly: false,
        },
        restore: {},
        saveAsImage: {
          name: "wgscovplot_variant_matrix",
        },
      },
    },
    grid: {
      containLabel: true,
      top: "5%"
    }
  };
}
