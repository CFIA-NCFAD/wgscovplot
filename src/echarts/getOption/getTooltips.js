import {meanCoverage, genomeCoverage, medianCoverage} from "../../coverageStat";
import {toTableHtml, getVariantComparison} from "../../util";

/**
 * Define options for tooltips
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {string} refSeq - Reference seq
 * @param {string} triggerOnType - mousemove or click
 * @param {boolean} isVariantSitesOnly - whether to show tooltips for variant sites only
 * @param {boolean} isVariantComparison - whether to compare variants across samples
 * @returns {Array<Object>}
 */
function getTooltips(samples, depths, variants, refSeq,
                     triggerOnType="mousemove", isVariantSitesOnly= false,
                     isVariantComparison= false) {
    return [
        {
            trigger: "axis",
            enterable: true,
            triggerOn: triggerOnType,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            position:"", // tooltip follows cursor
            formatter: function (params) {
                var output = "";
                var param = params[0];
                var i = param.axisIndex;
                if (i > samples.length) {
                    return output;
                }
                var sample = samples[i];
                var position = param.data[1];
                var depth = param.data[0];
                var zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
                var zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
                var meanCov = meanCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                var medianCov = medianCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                var genomeCov = genomeCoverage(depths, zoomStart, zoomEnd, i, 10).toFixed(2);
                var rows = [];
                const isVariantBar = params.find(element => {
                    if (element.componentSubType === "bar") {
                        return true;
                    }
                    return false
                });
                if (isVariantBar) {
                    if (isVariantComparison){
                        rows = getVariantComparison(samples, variants, depths, position, sample);
                    }
                    else {
                        rows = [
                            ["Position", position.toLocaleString()],
                            ["Coverage Depth", depth.toLocaleString()],
                        ];
                        Object.values(variants[i]).forEach(values => {
                            if (values['POS'] === position) {
                                for (const [key, value] of Object.entries(values)) {
                                    if (key !== 'POS' && key !== 'sample') {
                                        rows.push(
                                            ...[[key, value]]
                                        )
                                    }
                                }
                            }
                        })
                    }
                } else {
                    if (isVariantSitesOnly){
                        rows = [];
                    }
                    else{
                        rows = [
                            ["Position", position.toLocaleString()],
                            ["Coverage Depth", depth.toLocaleString()],
                        ];
                        rows.push(["Sequence", refSeq[position - 1]]);
                    }
                }
                if (rows.length){
                    output += "<h5>" + sample + "</h5>";
                    output += toTableHtml(["Position Info", ""], rows, "table small");
                    rows = [
                        [
                            "Range",
                            zoomStart.toLocaleString() + " - " + zoomEnd.toLocaleString(),
                        ],
                        ["Mean Coverage", meanCov + "X"],
                        ["Median Coverage", medianCov + "X"],
                        ["Genome Coverage ( >= 10x)", genomeCov + "%"],
                    ];
                    output += toTableHtml(["Coverage View Stats", ""], rows, "table small");
                }
                return output;
            },
        },
    ];
}

export {getTooltips};