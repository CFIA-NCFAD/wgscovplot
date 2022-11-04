import {toTableHtml} from "../../util";
import {getSegmentCoords, whichSegment} from "./getFluSegmentsInfo";
import {find, flatMap} from "lodash/collection";
import {uniqBy} from "lodash/array";
import {keys} from "lodash/object";
import {isEmpty} from "lodash/lang";
import {genomeCoverage, meanCoverage, medianCoverage} from "../../stats";
import {isNil} from "lodash";

/**
 * Get Variant sites comparison among samples
 * @param {string} sample - Sample of interest
 * @param {Array<string>} selectedSamples - Names of other samples being shown
 * @param {string} segment - Name of segment of interest
 * @param {Object} variants - Variant calling information in map of {sample: {segment: [variantCalls]}}
 * @param {Object} depths - Object of depths array
 * @param {Object} ref_seqs - The object of reference sequences of each segment
 * @param {Object} ref_ids - The object of reference ID of each segment
 * @param {number} position - Position in xAxis
 * @returns {Array<Array<string>>}
 */
function getFluVariantComparison(
    {
        sample,
        position,
        db: {
            selectedSamples,
            segment,
            depths,
            variants,
            ref_seqs,
            ref_ids,
        },
    }) {
    let variantArr = [];
    selectedSamples.forEach(s => {
        let variantCall = find(variants[s][segment],
            {"POS": position});
        if (!isNil(variantCall)) {
            variantArr.push(variantCall);
        } else {
            variantArr.push({
                "Sample": s,
                "POS": position,
                "REF_ID": ref_ids[s][segment],
                "Segment": segment,
                "Segment Length": depths[s][segment].length,
                "REF_SEQ": ref_seqs[s][segment][position - 1]
            });
        }
    });
    let unionKeys = uniqBy(flatMap(variantArr, keys));
    unionKeys.push("Coverage Depth"); // Add Coverage Depth row
    let rows = [];
    unionKeys.forEach(key => {
        let row = [];
        row.push(key);
        if (key === "Coverage Depth") {
            selectedSamples.forEach(s => {
                const segDepths = depths[s][segment];
                let rowValue;
                if (segDepths.length === 0) {
                    rowValue = `No result reported for segment ${segment}`;
                } else {
                    if (position <= segDepths.length) {
                        rowValue = segDepths[position - 1].toLocaleString();
                    } else {
                        //Out of range when segment length < padding array
                        rowValue = `No sequence at this position. Reference sequence 
                            ${ref_ids[s][segment]} is only 
                            ${segDepths.length} bp`;
                    }
                }
                row.push(rowValue);
            });
        } else {
            variantArr.forEach(element => {
                let variant = element[key];
                if (!isNil(variant)) {
                    if (key === "Sample" && variant === sample) {
                        // Bold highlight selected sample
                        row.push(variant.bold());
                    } else {
                        row.push(variant);
                    }
                } else {
                    // No information
                    row.push("");
                }
            });
        }
        rows.push(...[row]);
    });
    return rows;
}

/**
 * Get Coverage stat comparison among samples
 * @param {string[]} samples - Sample names
 * @param {string} segment - Segment name
 * @param {Object} depths - Object of depths array
 * @param {number} position - Position in xAxis
 * @param {number} low - low coverage threshold
 * @returns {(string[])[]}
 */
function getFluCoverageStatComparison(
    {
        samples,
        segment,
        depths,
        position,
        low
    }) {
    let rows = [];
    let tableHeader = [
        "Sample",
        `Depth at position ${position}`,
        "Segment",
        "Range",
        "Segment Length",
        "Mean Coverage (X)",
        "Median Coverage (X)",
        `Genome Coverage (>=${low}X) (%)`
    ];
    rows.push(...[tableHeader]);
    samples.forEach((sample) => {
        let depthArr = depths[sample][segment];
        let meanCov = meanCoverage(depthArr, 1, depthArr.length).toFixed(2);
        let medianCov = medianCoverage(depthArr, 1, depthArr.length).toFixed(2);
        let genomeCov = genomeCoverage(depthArr, 1, depthArr.length, low).toFixed(2);
        let coverageDepth;
        if (position <= depthArr.length) {
            coverageDepth = depthArr[position - 1].toLocaleString();
        } else if (depthArr.length === 0) {
            coverageDepth = `No result reported for segment ${segment}`;
        } else {
            coverageDepth = "No sequence at this position";
        }
        let row = [
            sample,
            coverageDepth,
            segment,
            `1-${depthArr.length}`,
            depthArr.length,
            meanCov,
            medianCov,
            genomeCov
        ];
        rows.push(row);
    });
    return rows;
}

/**
 * Get tooltips options
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Object} variants - The object of variants data
 * @param {Object} refSeqs - The object of reference sequences of each segment
 * @param {Object} refIDs - The object of reference ID of each segment
 * @param {string} triggerOnType - mousemove or click (tooltips options)
 * @param {boolean} showTooltipCovStats - whether to show coverage statistics (tooltips options)
 * @param {boolean} infoComparison - whether to compare variants/coverage stat across samples (tooltips options)
 * @param {number} low_coverage_threshold - low coverage threshold
 * @param {Object} primerData - The object of primer alignment of each segment
 * @returns {Array<Object>} - tooltips Option
 */
function tooltips(
    {
        samples,
        segments,
        depths,
        variants,
        refSeqs,
        refIDs,
        triggerOnType,
        infoComparison,
        showTooltipCovStats,
        low_coverage_threshold,
        primerData,
    }) {
    let segCoords = getSegmentCoords({samples, segments, depths});
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
                type: "line"
            },
            formatter: function (params) {
                let param = params[0];
                let output = "";
                // pos of xAxis in full scale
                let position = param.axisValue;
                let i = param.axisIndex;
                let positionRows = [];
                let coverageStatRows = [];
                let primerInfoRows = [];
                let sample = samples[i];
                // find segment index which pos belongs to
                let segment = whichSegment({position, segCoords});
                // get segment length
                let sequence = refSeqs[sample][segment];
                let segmentLength = sequence.length;
                // convert to pos in segment
                position = position - segCoords[segment].start + 1;
                let coverageDepth;
                let sampleSegDepths = depths[sample][segment];
                let refID = refIDs[sample][segment];
                if (segmentLength === 0) {
                    coverageDepth = `No result reported for segment ${segment}`;
                } else {
                    if (position <= segmentLength) {
                        // get coverage depth for pos
                        coverageDepth = sampleSegDepths[position - 1].toLocaleString();
                    } else {
                        coverageDepth = `No sequence at this position. 
                        Reference sequence ${refID} 
                        is only ${sampleSegDepths.length} bp`;
                    }
                }
                // check if current position is variant sites
                const variantBar = params.find(element => element.componentSubType === "bar");
                if (variantBar) { //tooltips at Variant sites
                    if (infoComparison) {
                        positionRows = getFluVariantComparison(
                            {
                                sample,
                                samples,
                                segment,
                                depths,
                                variants,
                                refSeqs,
                                refIDs,
                                position,
                            });
                    } else {
                        let foundObj = find(variants[samples[param.axisIndex]][segment],
                            {"POS": position});
                        if (foundObj !== undefined && foundObj !== null) {
                            for (const [key, value] of Object.entries(foundObj)) {
                                if (key !== "Sample") { // do not write row for sample name
                                    positionRows.push(
                                        ...[[key, value]]
                                    );
                                }
                            }
                        }
                        positionRows.push(["Coverage Depth", coverageDepth]);
                    }
                } else { // tooltips for Non-Variant Sites
                    if (position > segmentLength) { //Out of range when segment length < padding array
                        positionRows = [
                            ["Position", position],
                            [coverageDepth, ""]
                        ];
                        output += `<h5>Sample: ${sample}</h5>`;
                        output += toTableHtml(
                            {
                                headers: ["Position Info", ""],
                                rows: positionRows,
                                classes: "table small"
                            });
                        return output;
                    } // Pos within segment length
                    positionRows = [
                        ["Segment", segment],
                        ["Segment Length", segmentLength],
                        ["REF_ID", refID],
                        ["POS", position],
                        ["REF_SEQ", sequence[position - 1]],
                        ["ALT_SEQ", ""],
                        ["ALT_FREQ", ""],
                        ["Coverage Depth", coverageDepth]
                    ];
                }
                if (positionRows.length) { // write rows to table
                    output += `<h5>Sample: ${sample}</h5>`;
                    output += toTableHtml(
                        {
                            headers: ["Position Info", ""],
                            rows: positionRows,
                            classes: "table small"
                        });
                    if (showTooltipCovStats) {
                        if (infoComparison) {
                            coverageStatRows = getFluCoverageStatComparison(
                                {
                                    samples,
                                    segment,
                                    depths,
                                    position,
                                    low: low_coverage_threshold,
                                });
                        } else {
                            let meanCov = meanCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
                            let medianCov = medianCoverage(sampleSegDepths, 1, segmentLength).toFixed(2);
                            let genomeCov = genomeCoverage(sampleSegDepths, 1, segmentLength, low_coverage_threshold).toFixed(2);
                            coverageStatRows = [
                                ["Mean Coverage", `${meanCov}X`],
                                ["Median Coverage", `${medianCov}X`],
                                [`Genome Coverage (>= ${low_coverage_threshold}X)`, `${genomeCov}%`],
                            ];
                        }
                        output += toTableHtml(
                            {
                                headers: ["Coverage View Stats", ""],
                                rows: coverageStatRows,
                                classes: "table small"
                            });
                    }
                    if (!isEmpty(primerData)) {
                        let primerInfo = primerData[sample][segment];
                        for (let i = 0; i < primerInfo.length; i++) {
                            let {
                                start,
                                end,
                                name,
                                query_aligned,
                                matched_aligned,
                                target_aligned,
                                edit_distance,
                                cigar,
                                other_locations,
                            } = primerInfo[i];
                            if (position >= start && position <= end) {
                                primerInfoRows = [
                                    ["Primer Name", name],
                                    ["Primer Sequence", query_aligned],
                                    ["Match Aligned", matched_aligned],
                                    ["Ref Sequence", target_aligned],
                                    ["Cigar", cigar],
                                    ["Start", start + 1],
                                    ["End", end + 1],
                                    ["Edit Distance", edit_distance],
                                    ["Other Locations", other_locations]
                                ];
                            }
                        }
                        if (primerInfoRows.length) {
                            output += toTableHtml({
                                headers: ["Primer Info", ""],
                                rows: primerInfoRows,
                                classes: "table small"
                            });
                        }
                    }
                    return output;
                }
            }
        }
    ];
}

export {tooltips};