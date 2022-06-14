import {toTableHtml} from "../../util";
import {map, filter, orderBy, find} from "lodash/collection";
import {union, uniq, reverse} from "lodash/array";

/**
 * Get tooltip for variant heatmap
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} mutations - Array of mutation name
 * @param {Array<Array<Object>>} variants - The object of variants data
 * @returns Array<Object>
 */
function getTooltipHeatmap(samples, mutations, variants) {

    let toolTips = [
        {
            enterable: true,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            formatter: function (params) {
                let output = '';
                let mutationName = '';
                let sample = '';
                // values may be undefined when mouse move and zoom/in out heat map (causing error log)
                if (params.value !== undefined && params.value !== null) {
                    mutationName = mutations[params.value[0]];
                    sample = samples[params.value[1]];
                }
                let rows = [];
                if (variants[sample] !== undefined && variants[sample] !== null) {
                    let foundObj = find(Object.values(variants[sample]), {"mutation": mutationName});
                    if (foundObj !== undefined && foundObj !== null) {
                        for (const [key, value] of Object.entries(foundObj)) {
                            rows.push(...[[key, value]]);
                        }
                    }
                    if (rows.length) {
                        output += toTableHtml(["", ""], rows, "table small");
                    } else {
                        output += ("Not detected mutation " + mutationName.bold() + " in sample " + sample.bold());
                    }
                } else {
                    output += "Sample " + sample.bold() + " has no variant information";
                }
                return output;
            }
        }
    ];
    return toolTips;
}

/**
 * Prepare data for Variant heatmap
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<Object>>} variants - The object of variants data
 * @returns {Array<>} - Array of mutation name and array of alt frequency matrix
 */
function getMutationMatrix(samples, variants) {
    let samplesObject = [];
    samples.forEach(sample => {
        let sampleInfo = variants[sample];
        samplesObject = union(samplesObject, filter(sampleInfo, "mutation"));
    });
    let samplesObjectSorted = orderBy(samplesObject, "POS", "asc");
    let mutation = uniq(map(samplesObjectSorted, "mutation"));
    let altFreqMatrix = [];
    for (let i = 0; i < samples.length; i++){
         for (let j = 0;  j < mutation.length; j++){
              let foundObj = find(samplesObjectSorted, {"sample": samples[i], "mutation": mutation[j]});
              if (foundObj !== undefined && foundObj !== null){
                 altFreqMatrix.push([j, samples.length - 1 - i, foundObj.ALT_FREQ]);
              }else{
                 altFreqMatrix.push([j, samples.length - 1 - i, 0]);
              }
         }
    }
    return [mutation, altFreqMatrix];
}
/**
 * Define all options for variant heatmap
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} mutations - Array of mutation name
 * @param {Array<Array<number>>} variantMatrix - Array of vairant heatmap value
 * @param {Array<Array<Object>>} variants - The object of variants data
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
 3: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G4181T(orf1ab:A1306S)', POS: 4181, REF: 'G', …}
 4: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A6022G(orf1ab:P1919P)', POS: 6022, REF: 'A', …}
 5: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C6286T(orf1ab:T2007T)', POS: 6286, REF: 'C', …}
 6: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C6402T(orf1ab:P2046L)', POS: 6402, REF: 'C', …}
 7: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C7124T(orf1ab:P2287S)', POS: 7124, REF: 'C', …}
 8: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C8956T(orf1ab:Y2897Y)', POS: 8956, REF: 'C', …}
 9: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C8986T(orf1ab:D2907D)', POS: 8986, REF: 'C', …}
 10: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G9053T(orf1ab:V2930L)', POS: 9053, REF: 'G', …}
 11: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C10029T(orf1ab:T3255I)', POS: 10029, REF: 'C', …}
 12: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A11201G(orf1ab:T3646A)', POS: 11201, REF: 'A', …}
 13: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A11332G(orf1ab:V3689V)', POS: 11332, REF: 'A', …}
 14: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G13897T(orf1ab:V4544V)', POS: 13897, REF: 'G', …}
 15: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C14408T(orf1ab:L4715L)', POS: 14408, REF: 'C', …}
 16: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G15451A(orf1ab:A5062A)', POS: 15451, REF: 'G', …}
 17: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C15848T(orf1ab:L5195L)', POS: 15848, REF: 'C', …}
 18: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C16466T(orf1ab:H5401Y)', POS: 16466, REF: 'C', …}
 19: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A17236G(orf1ab:L5657L)', POS: 17236, REF: 'A', …}
 20: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C18657T(orf1ab:P6131L)', POS: 18657, REF: 'C', …}
 21: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C19220T(orf1ab:L6319L)', POS: 19220, REF: 'C', …}
 22: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C21618G(S:T19R)', POS: 21618, REF: 'C', …}
 23: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G21987A(S:G142D)', POS: 21987, REF: 'G', …}
 24: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'GAGTTCA22028G(S:E156_R158DELINSG)', POS: 22028, REF: 'GAGTTCA', …}
 25: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'T22917G(S:L452R)', POS: 22917, REF: 'T', …}
 26: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C22995A(S:T478K)', POS: 22995, REF: 'C', …}
 27: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A23403G(S:D614G)', POS: 23403, REF: 'A', …}
 28: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C23604G(S:P681R)', POS: 23604, REF: 'C', …}
 29: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C24208T(S:I882I)', POS: 24208, REF: 'C', …}
 30: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A24253T(S:P897P)', POS: 24253, REF: 'A', …}
 31: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G24410A(S:D950N)', POS: 24410, REF: 'G', …}
 32: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C25469T(ORF3a:S26L)', POS: 25469, REF: 'C', …}
 33: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C26256T(E:F4F)', POS: 26256, REF: 'C', …}
 34: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'T26767C(M:I82T)', POS: 26767, REF: 'T', …}
 35: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C27005T(M:I161I)', POS: 27005, REF: 'C', …}
 36: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'T27638C(ORF7a:V82A)', POS: 27638, REF: 'T', …}
 37: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C27752T(ORF7a:T120I)', POS: 27752, REF: 'C', …}
 38: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C27874T', POS: 27874, REF: 'C', …}
 39: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C27999T(ORF8:P36S)', POS: 27999, REF: 'C', …}
 40: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'AGATTTC28247A(ORF8:D119_F120del)', POS: 28247, REF: 'AGATTTC', …}
 41: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'C28253A(ORF8:F120L)', POS: 28253, REF: 'C', …}
 42: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'TA28270T', POS: 28270, REF: 'TA', …}
 43: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'A28461G(N:D63G)', POS: 28461, REF: 'A', …}
 44: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G28881T(N:R203M)', POS: 28881, REF: 'G', …}
 45: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G28916T(N:G215C)', POS: 28916, REF: 'G', …}
 46: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G29402T(N:D377Y)', POS: 29402, REF: 'G', …}
 47: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G29742T', POS: 29742, REF: 'G', …}
 48: {sample: 'SRR17230671', CHROM: 'MN908947.3', mutation: 'G29747T', POS: 29747, REF: 'G', …}
 SRR17230672: (54) [{…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}]
 SRR17230675: (56) [{…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}]
 SRR17230678: (50) [{…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}]
 SRR17230680: (43) [{…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}, {…}
 */

function getVariantHeatmapOption(samples, variants) {
    let mutationMatrixInfo = getMutationMatrix(samples, variants);
    let dataSamples = [...samples].reverse();
    let chartOptions = {
        xAxis: {
            type: "category",
            data: mutationMatrixInfo[0],
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
                show: true
            }
        ],
        visualMap: {
            min: 0.00,
            max: 1.00,
            precision: 2,
            calculable: true,
            orient: 'vertical',
            left: 'right',
            top: '5%',
            inRange: {
                color: [
                    '#a50026',
                    '#f46d43',
                    '#fdae61',
                    '#fee08b',
                    '#f7ff00',
                    '#006837'
                ]
            }
        },
        series: [
            {
                type: "heatmap",
                data: mutationMatrixInfo[1],
                label: {
                    show: false
                },
                emphasis: {
                    itemStyle: {
                        shadowBlur: 10,
                        shadowColor: 'rgba(0, 0, 0, 0.5)'
                    }
                }
            }
        ],
        tooltip: getTooltipHeatmap(dataSamples, mutationMatrixInfo[0], variants),
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
    return chartOptions;
}

export {getVariantHeatmapOption};