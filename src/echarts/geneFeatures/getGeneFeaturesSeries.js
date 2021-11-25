import {renderGeneFeatures} from "./renderGeneFeatures";
import {toTableHtml} from "../../util";

/**
 * Define options for gene features charts
 * @param {number} index - gene feature is displayed in the last index of grid
 * @returns {Array<Dict[]>}
 */
function getGeneFeatureSeries(index) {
    var featureSeries = [];
    if (amplicon === 'True' || geneFeature === 'True'){
        featureSeries.push({
            type: "custom",
            xAxisIndex: index,
            yAxisIndex: index,
            renderItem: renderGeneFeatures,
            labelLayout: {
                hideOverlap: false,
            },
            data: geneFeatureData,
            tooltip: {
                trigger: "item",
                enterable: true,
                appendToBody: true,
                renderMode: "html",
                borderRadius: 6,
                borderWidth: 2,
                showContent: "true",
                position: 'top',
                textStyle: {
                    fontSize: 15,
                    fontWeight: "bolder",
                },
                formatter: function (params) {
                    var output = "";
                    var rows = [
                        [
                            "Range",
                            params.value.start.toLocaleString() + " - " + params.value.end.toLocaleString(),
                        ],
                        ["Length", (params.value.end - params.value.start + 1).toLocaleString()],
                        ["Strand", params.value.strand],
                    ];
                    output += toTableHtml([params.name, ""], rows, "table small");
                    return output;
                },
            },
        });
    }
    return featureSeries;
}

export {getGeneFeatureSeries};