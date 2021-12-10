import {meanCoverage, genomeCoverage, medianCoverage} from "../../coverageStat";
import {toTableHtml, getVariantComparison, getCoverageStatComparison} from "../../util";


/**
 * Define options for tooltips
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {string} refSeq - Reference seq
 * @param {string} triggerOnType - mousemove or click
 * @returns {Array<Object>}
 */
function getTooltips(samples, depths, variants, refSeq,
                     triggerOnType="mousemove", isVariantComparison = false) {
    return [
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
                var output = "";
                var param = params[0];
                var i = param.axisIndex;
                if (i > samples.length) {
                    return output;
                }
                var sample = samples[i];
                var position = param.axisValue;
                var depth = depths[i][position-1];
                var zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
                var zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
                var positionRows = [];
                var coverageStatRows = [];
                const isVariantBar = params.find(element => {
                    if (element.componentSubType === "bar") {
                        return true;
                    }
                    return false
                });
                if (isVariantBar) {
                    if (isVariantComparison){
                        positionRows = getVariantComparison(samples, variants, depths, position, sample);
                    }
                    else {
                        positionRows = [
                            ["Position", position.toLocaleString()],
                            ["Coverage Depth", depth.toLocaleString()],
                        ];
                        Object.values(variants[i]).forEach(values => {
                            if (values['POS'] === position) {
                                for (const [key, value] of Object.entries(values)) {
                                    if (key !== 'POS' && key !== 'sample') {
                                        positionRows.push(
                                            ...[[key, value]]
                                        )
                                    }
                                }
                            }
                        })
                    }
                } else {
                    positionRows = [
                        ["Position", position.toLocaleString()],
                        ["Coverage Depth", depth.toLocaleString()],
                    ];
                    positionRows.push(["Sequence", refSeq[position-1]]);
                }
                if (positionRows.length){
                    output += "<h5>" + "Selected sample: " + sample + "</h5>";
                    output += toTableHtml(["Position Info", ""], positionRows, "table table-hover table-bordered table-responsive-md");
                    if (isVariantComparison){
                        coverageStatRows = getCoverageStatComparison(samples, depths, zoomStart, zoomEnd, sample);
                    }
                    else{
                        var meanCov = meanCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                        var medianCov = medianCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                        var genomeCov = genomeCoverage(depths, zoomStart, zoomEnd, i, 10).toFixed(2);
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
                    output += toTableHtml(["Coverage View Stats", ""], coverageStatRows, "table table-hover table-bordered table-responsive-md");
                }
                return output;
            },
        },
    ];
}

export {getTooltips};