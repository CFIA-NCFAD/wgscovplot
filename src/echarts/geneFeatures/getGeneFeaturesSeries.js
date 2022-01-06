import {getGeneFeatureRenderer} from "./getGeneFeatureRenderer";
import {toTableHtml} from "../../util";

/**
 * Define options for gene features charts
 * @param {Array<Object>} geneFeatureData - Array of dictionary geneFeature or amplicon data
 * @param {number} index - gene feature is displayed in the last index of grid
 * @param {boolean} geneFeature - whether to plot gene feature or not (true for false)
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @returns {Array<Object>}
 */
function getGeneFeatureSeries(geneFeatureAmpliconData, index, geneFeature, amplicon) {
    var featureSeries = [];
    if (amplicon || geneFeature){
        featureSeries.push({
            type: "custom",
            xAxisIndex: index,
            yAxisIndex: index,
            renderItem: getGeneFeatureRenderer(true, geneFeatureAmpliconData),
            labelLayout: {
                hideOverlap: false,
            },
            data: geneFeatureAmpliconData,
            tooltip: {
                trigger: "item",
                enterable: true,
                appendToBody: true,
                triggerOn: "mousemove",
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