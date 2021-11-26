import {ntColor} from "../../util";

/**
 * Define options for variant bar charts
 * @param {Array<Dict[]>} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {string} refSeq - Reference seq
 * @returns {Array<Dict[]>}
 */
function getVariantsSeries(variants, depths, refSeq) {
    var variantsSeries = [];
    for (var [i, varMap] of variants.entries()) {
        (function (i, varMap) {
            variantsSeries.push({
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

            });
        })(i, varMap);
    };
    return variantsSeries;
}

export {getVariantsSeries};