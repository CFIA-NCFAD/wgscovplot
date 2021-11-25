import {meanCoverage, genomeCoverage, medianCoverage} from "../../coverageStat";
import {toTableHtml} from "../../ulti";

/**
 * Define options for tooltips
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @returns {Array<Dict[]>}
 */
function getTooltips(samples, depths, variants) {
    return [
        {
            trigger: "axis",
            enterable: true,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            position: function (pos, params, dom, rect, size) {
                // tooltip will be fixed on the right if mouse hovering on the left,
                // and on the left if hovering on the right.
                var obj = {top: 5};
                obj[['left', 'right'][+(pos[0] < size.viewSize[0] / 2)]] = 5;
                return obj;
            },
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
                output += "<h5>" + sample + "</h5>";
                var rows = [
                    ["Position", position.toLocaleString()],
                    ["Depth", depth.toLocaleString()],
                ];
                if (params.length > 1) {
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
                } else {
                    rows.push(["Sequence", window.refSeq[position - 1]]);
                }
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

                return output;
            },
        },
    ];
}

export {getTooltips};