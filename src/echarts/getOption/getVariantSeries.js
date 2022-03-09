import {ntColor} from "../../util";

/**
 * Define options for variant bar charts
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {string} refSeq - Reference seq
 * @param {boolean} isVariantSites- whether to show tooltips for variant sites
 * @param {boolean} isShowMutation - whether to show Muation below Variant Sites
 * @returns {Array<Object>}
 */
function getVariantsSeries(variants, depths, refSeq, isVariantSites,
                           isShowMutation) {
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
                label: {
                    show: isShowMutation,
                    position: "bottom",
                    align: "left",
                    verticalAlign: "middle",
                    distance: 10,
                    color: "inherit",
                    rotate: -30,
                    formatter: function (params){
                        var pos = params.data[0];
                        var output = "";
                        Object.values(variants[i]).forEach(values => {
                            if (values['POS'] === pos) {
                                if (values["mutation"] !== undefined && values["mutation"] !== null){
                                    output += values["mutation"];
                                }
                            }
                        });
                        return output;
                    }
                },
                labelLayout:{
                    hideOverlap: true
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