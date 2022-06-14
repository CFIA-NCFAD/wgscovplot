import {meanCoverage, genomeCoverage, medianCoverage} from "../../coverageStat";
import {toTableHtml, getVariantComparison, getCoverageStatComparison} from "../../util";
import {find} from "lodash/collection";


/**
 * Define options for tooltips
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {string} refSeq - Reference seq
 * @param {string} triggerOnType - mousemove or click
 * @param {boolean} isInfoComparison - whether to compare variants/ Coverage Stat across samples
 * @param {boolean} isCovergateStatView - whether to show Coverage stat
 * @returns {Array<Object>}
 */
function getTooltips(samples, depths, variants, refSeq,
                     triggerOnType="mousemove", isInfoComparison=true,
                     isCovergateStatView=false) {
    let toolTips = [
        {
            trigger: "axis",
            enterable: true,
            triggerOn: triggerOnType,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            position: "cursor",
            axisPointer: {
                type: 'line'
            },
            formatter: function (params) {
                let output = "";
                let param = params[0];
                let i = param.axisIndex;
                if (i > samples.length) {
                    return output;
                }
                let sample = samples[i];
                let position = param.axisValue;
                let depth = depths[i][position-1];
                let zoomStart = 1;
                let zoomEnd = refSeq.length;
                if (isCovergateStatView){
                    zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
                    zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
                }
                let positionRows = [];
                let coverageStatRows = [];
                const isVariantBar = params.find(element => {
                    if (element.componentSubType === "bar") {
                        return true;
                    }
                    return false;
                });
                if (isVariantBar) {
                    if (isInfoComparison){
                        positionRows = getVariantComparison(samples, variants, depths, position, sample);
                    }
                    else {
                        positionRows = [
                            ["Position", position.toLocaleString()],
                            ["Coverage Depth", depth.toLocaleString()],
                        ];
                        let foundObj = find(Object.values(variants[i]), {"POS": position});
                        if (foundObj !== undefined && foundObj !== null){
                            for (const [key, value] of Object.entries(foundObj)) {
                                if (key !== 'POS' && key !== 'sample') {
                                    positionRows.push(
                                        ...[[key, value]]
                                    );
                                }
                            }
                        }
                    }
                } else {
                    positionRows = [
                        ["Position", position.toLocaleString()],
                        ["Coverage Depth", depth.toLocaleString()],
                    ];
                    positionRows.push(["Sequence", refSeq[position-1]]);
                }
                if (positionRows.length){
                    output += "<h5>" + "Sample: " + sample + "</h5>";
                    output += toTableHtml(["Position Info", ""], positionRows, "table small");
                    if (isCovergateStatView){
                        if (isInfoComparison){
                            coverageStatRows = getCoverageStatComparison(samples, depths, zoomStart, zoomEnd, sample, position);
                        }
                        else{
                            let meanCov = meanCoverage(depths[i], zoomStart, zoomEnd).toFixed(2);
                            let medianCov = medianCoverage(depths[i], zoomStart, zoomEnd).toFixed(2);
                            let genomeCov = genomeCoverage(depths[i], zoomStart, zoomEnd, 10).toFixed(2);
                            coverageStatRows = [
                                [
                                    "Range",
                                    zoomStart.toLocaleString() + " - " + zoomEnd.toLocaleString(),
                                ],
                                ["Mean Coverage", meanCov + "X"],
                                ["Median Coverage", medianCov + "X"],
                                ["Genome Coverage ( >= 10x)", genomeCov + "%"],
                            ];
                        }
                        output += toTableHtml(["Coverage View Stats", ""], coverageStatRows, "table small");
                    }
                }
                return output;
            },
        },
    ];
    return toolTips;
}

export {getTooltips};