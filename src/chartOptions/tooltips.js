import {find, flatMap} from "lodash/collection";
import {get, isNil} from "lodash";
import {uniqBy} from "lodash/array";
import {keys} from "lodash/object";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {toFloat32Array, toTableHtml} from "../util";

/**
 * Function get Variant Comparison across samples
 * @param {WgsCovPlotDB} db
 * @param {number} position - x-axis position
 * @param {string} sampleInFocus - sample
 * @returns <Array<Array<string>> - Variant comparison across samples
 */
function getVariantComparison(
    {
        db,
        position,
        sampleInFocus = "",
    }) {
    const {
        selectedSamples,
        variants,
        depths,
    } = db;
    let rows = [];
    let variantArr = [];
    for (let sample of selectedSamples) {
        let sampleVariants = variants[sample];
        let variant = find(sampleVariants, {POS: position + ""})
        if (!isNil(variant)) {
            variantArr.push(variant);
        } else {
            variantArr.push({sample, POS: position});
        }
    }
    let unionKeys = uniqBy(flatMap(variantArr, keys));
    unionKeys.push("Coverage Depth"); // Add Coverage Depth row
    unionKeys.forEach(key => {
        let row = [];
        row.push(key);
        if (key === "Coverage Depth") {
            for (let sample of selectedSamples) {
                row.push(toFloat32Array(depths[sample])[position - 1].toLocaleString());
            }
        } else {
            for (let variant of variantArr) {
                let val = get(variant, key, "");
                row.push((val === sampleInFocus) ? val.bold() : val);
            }
        }
        rows.push(...[row]);
    });
    return rows;
}


/**
 * Function get Coverage Stat comparison across samples
 * @param {WgsCovPlotDB} db
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} low_coverage_threshold - low coverage threshold
 * @param {number} position - Selected position
 * @returns <Array<Array<string>> - Coverage Stat comparison across samples
 */
function getCoverageStatComparison(
    {
        db,
        start,
        end,
        position
    }) {
    const {
        selectedSamples,
        depths,
        low_coverage_threshold,
    } = db;
    let rows = [];
    let tableHeader = [
        "Sample",
        `Depth at position ${position.toLocaleString()}`,
        "Range",
        "Mean Coverage (X)",
        "Median Coverage (X)",
        `Genome Coverage (>=${low_coverage_threshold}X) (%)`,
    ];
    rows.push(...[tableHeader]);
    for (let sample of selectedSamples) {
        let sampleDepths = toFloat32Array(depths[sample]);
        let meanCov = meanCoverage(sampleDepths, start, end).toFixed(2);
        let medianCov = medianCoverage(sampleDepths, start, end).toFixed(2);
        let genomeCov = genomeCoverage(sampleDepths, start, end, low_coverage_threshold).toFixed(2);
        let coverageDepth = sampleDepths[position - 1];
        let row = [
            sample,
            coverageDepth.toLocaleString(),
            `${start.toLocaleString()}-${end.toLocaleString()}`,
            meanCov,
            medianCov,
            genomeCov,
        ];
        rows.push(...[row]);
    }
    return rows;
}

function tooltipFormatter({db}) {
    return function (params) {
        const {
            selectedSamples,
            depths,
            chart,
            crossSampleComparisonInTooltips,
            variants,
            ref_seq,
            low_coverage_threshold,
            showCovStatsInTooltips,
        } = db;
        let output = "";
        let [{
            axisIndex,
            axisValue: position
        }] = params;
        if (axisIndex >= selectedSamples.length) {
            return output;
        }
        let sample = selectedSamples[axisIndex];
        let sampleDepths = toFloat32Array(depths[sample]);
        let depth = sampleDepths[position - 1];
        //let depth1 = window.atob(depths[sample]).charCodeAt(position - 1).buffer
        //console.log(depth1)
        let [dataZoom] = chart.getOption().dataZoom;
        let zoomStart = Math.floor(dataZoom.startValue);
        let zoomEnd = Math.floor(dataZoom.endValue);
        let positionRows = [];
        let coverageStatRows = [];
        const isVariantBar = params.find(element => {
            return element.componentSubType === "bar";

        });
        if (isVariantBar) {
            if (crossSampleComparisonInTooltips) {
                positionRows = getVariantComparison({
                    db,
                    position,
                    sampleInFocus: sample,
                });
            } else {
                positionRows = [
                    ["Position", position.toLocaleString()],
                    ["Coverage Depth", depth.toLocaleString()],
                ];
                let foundObj = find(variants[sample], {"POS": position.toString()}, 0);
                if (!isNil(foundObj)) {
                    for (const [key, value] of Object.entries(foundObj)) {
                        if (key !== "POS" && key !== "sample") {
                            positionRows.push([key, value]);
                        }
                    }
                }
            }
        } else {
            positionRows = [
                ["Pos", position.toLocaleString()],
                ["Coverage Depth", depth.toLocaleString()],
            ];
            positionRows.push(["Sequence", ref_seq[position - 1]]);
        }
        if (positionRows.length) {
            output += `<h5>Sample: ${sample}</h5>`;
            output += toTableHtml({
                headers: ["Position Info", ""],
                rows: positionRows,
            });
            if (showCovStatsInTooltips) {
                if (crossSampleComparisonInTooltips) {
                    coverageStatRows = getCoverageStatComparison({
                        db,
                        start: zoomStart,
                        end: zoomEnd,
                        position
                    });
                } else {
                    let meanCov = meanCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
                    let medianCov = medianCoverage(sampleDepths, zoomStart, zoomEnd).toFixed(2);
                    let genomeCov = genomeCoverage(sampleDepths, zoomStart, zoomEnd, low_coverage_threshold).toFixed(2);
                    coverageStatRows = [
                        [
                            "Range",
                            `${zoomStart.toLocaleString()} - ${zoomEnd.toLocaleString()}`,
                        ],
                        ["Mean Coverage", `${meanCov}X`],
                        ["Median Coverage", `${medianCov}X`],
                        [`Genome Coverage (>= ${low_coverage_threshold}X)`, `${genomeCov}%`],
                    ];
                }
                output += toTableHtml(
                    {
                        headers: ["Coverage View Stats", ""],
                        rows: coverageStatRows,
                    });
            }
        }
        return output;
    };
}

/**
 * Define options for tooltips
 * @param {Object} db - wgscovplot DB object
 * @returns {Array<Object>}
 */
function getTooltips(db) {
    const {
        tooltipTriggerOn = "mousemove",
    } = db;
    return [
        {
            trigger: "axis",
            enterable: true,
            triggerOn: tooltipTriggerOn,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            position: "cursor",
            axisPointer: {
                type: "line"
            },
            formatter: tooltipFormatter({db}),
        },
    ];
}

export {
    getTooltips
};
