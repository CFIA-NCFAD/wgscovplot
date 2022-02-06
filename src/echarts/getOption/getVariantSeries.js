import {ntColor} from "../../util";

/**
 * Define options for variant bar charts
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {boolean} isVariantSites- whether to show tooltips for variant sites
 * @param {string} refSeq - Reference seq
 * @returns {Array<Object>}
 */
function getVariantsSeries(variants, depths, refSeq, isVariantSites) {
    var variantSeries = [];
    for (var [i, varMap] of variants.entries()) {
        (function (i, varMap) {
            variantSeries.push({
                type: "bar",
                xAxisIndex: i,
                yAxisIndex: i,
                data: Object.values(varMap).map((x) => [parseInt(x['POS']), depths[i][x['POS']]]),
                barWidth: 2,
                itemStyle: {
                    color: function (params) {
                        var pos = params.data[0];
                        var nt = refSeq[pos-1]
                        if (ntColor.hasOwnProperty(nt)) {
                            return ntColor[nt];
                        }
                        return "#333";
                    },
                },
                tooltip:{
                    trigger: isVariantSites ? "axis" : "none"
                }
            });
        })(i, varMap);
    }
    return variantSeries;
}


export {getVariantsSeries};