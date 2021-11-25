
/**
 * Dict of color for coloring variant position in the chart
 * @type {Dict[]}
 */
const ntColor = {
    "A": "#ea5e48",
    "C": "#eaca48",
    "G": "#6ad82b",
    "T": "#2b87d8",
}

/**
 * Define options for variant bar charts
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @returns {Array<Dict[]>}
 */
function getVariantsSeries(variants, depths) {
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
                        var nt = window.refSeq[pos-1]
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