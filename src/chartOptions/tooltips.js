import {find, flatMap} from "lodash/collection";
import {get, isNil} from "lodash";
import {uniqBy} from "lodash/array";
import {keys} from "lodash/object";
import {genomeCoverage, meanCoverage, medianCoverage} from "../stats";
import {toFloat32Array} from "../util";

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

export {
    getVariantComparison,
    getCoverageStatComparison
};
