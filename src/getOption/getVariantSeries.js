import {ntColor} from "../util";

/**
 * Define options for variant bar charts
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {string} refSeq - Reference seq
 * @param {boolean} variantSites- whether to show tooltips for variant sites
 * @param {boolean} showMutation - whether to show Mutation below Variant Sites
 * @param {boolean} hideOverlapMutation - whether to hide overlapping mutation under variants sites
 * @returns {Array<Object>}
 */
function getVariantsSeries(variants, depths, refSeq, variantSites,
                           showMutation, hideOverlapMutation) {
    let variantSeries = [];
    for (let [i, varMap] of variants.entries()) {
        (function (i, varMap) {
            variantSeries.push({
                type: "bar",
                xAxisIndex: i,
                yAxisIndex: i,
                data: Object.values(varMap).map((x) => [parseInt(x.POS), depths[i][parseInt(x.POS) - 1]]),
                barWidth: 2,
                itemStyle: {
                    color: function (params) {
                        let pos = params.data[0];
                        let nt = refSeq[pos - 1];
                        if (Object.prototype.hasOwnProperty.call(ntColor, nt)) {
                            return ntColor[nt];
                        }
                        return "#333";
                    },
                },
                label: {
                    show: showMutation,
                    position: "bottom",
                    align: "left",
                    verticalAlign: "middle",
                    distance: 10,
                    color: "inherit",
                    rotate: -30,
                    formatter: function (params) {
                        let pos = params.data[0];
                        let output = "";
                        Object.values(variants[i]).forEach(values => {
                            if (values.POS === pos) {
                                output += (values.REF + values.POS + values.ALT);
                            }
                        });
                        return output;
                    }
                },
                labelLayout: {
                    hideOverlap: hideOverlapMutation
                },
                tooltip: {
                    trigger: variantSites ? "axis" : "none"
                }
            });
        })(i, varMap);
    }
    return variantSeries;
}


export {getVariantsSeries};