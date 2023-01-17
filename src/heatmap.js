import {toTableHtml} from "./util";
import {filter, find, isEmpty, isNil, map, orderBy, union, uniq} from "lodash";

/**
 * Get tooltip for variant heatmap
 * @param {Array<string>} selectedSamples - An array of samples name
 * @param {Array<string>} mutations - Array of mutation name
 * @param {SampleVariantCalls} variants - The object of variants data
 * @returns {Object}
 */
function getTooltipHeatmap(
    {
        dataSamples,
        mutations,
        variants,
    }) {
    return {
        enterable: true,
        appendToBody: true,
        renderMode: "html",
        showContent: true,
        confine: true,
        formatter: function ({value}) {
            let output = "";
            let mutation = "";
            let sample = "";
            // values may be undefined when mouse move and zoom/in out heat map (causing error log)
            if (!isNil(value)) {
                mutation = mutations[value[0]];
                sample = dataSamples[value[1]];
            }
            let rows = [];
            let sampleVariants = variants[sample];
            if (!isNil(sampleVariants)) {
                let variantOfInterest = find(Object.values(sampleVariants), {mutation});
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

/**
 * Prepare data for Variant heatmap
 * @param {String[]} samples - An array of samples name
 * @param {SampleVariantCalls|SampleSegmentVariantCalls} variants - The object of variants data
 * @returns {[string[], number[]} - Tuple of mutation names array and array of mutation AF values for heatmap
 */
function getMutationMatrix(samples, variants) {
    let sampleVariants = [];
    samples.forEach(sample => {
        let variantCalls = variants[sample];
        sampleVariants = union(sampleVariants, filter(variantCalls, "mutation"));
    });
    sampleVariants = orderBy(sampleVariants, (obj) => parseInt(obj.POS));
    let mutations = uniq(map(sampleVariants, "mutation"));
    let altFreqMatrix = [];
    for (let i = 0; i < samples.length; i++) {
        let sample = samples[i];
        for (let j = 0; j < mutations.length; j++) {
            let mutation = mutations[j];
            let foundObj = find(sampleVariants, {"sample":sample, "mutation":mutation});
            altFreqMatrix.push([
                j,
                samples.length - 1 - i,
                !isNil(foundObj) ? parseFloat(foundObj.ALT_FREQ) : 0
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
function getVariantHeatmapOption(db) {
    const {
        selectedSamples,
        variants,
    } = db;
    let [mutations, matrix] = getMutationMatrix(selectedSamples, variants);
    let dataSamples = [...selectedSamples].reverse();
    let chartOption;
    chartOption = {
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
                show: true
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
                    "#a50026",
                    "#f46d43",
                    "#fdae61",
                    "#fee08b",
                    "#f7ff00",
                    "#006837"
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
        tooltip: getTooltipHeatmap({dataSamples, mutations, variants}),
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
    return chartOption;
}

export {getVariantHeatmapOption};