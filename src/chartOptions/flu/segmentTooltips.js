import {toTableHtml} from "../../util";
import {whichSegment} from "./segmentInfo";
import {find, flatMap} from "lodash/collection";
import {uniqBy} from "lodash/array";
import {keys} from "lodash/object";
import {genomeCoverage, meanCoverage, medianCoverage} from "../../stats";
import {has, isNil} from "lodash";

/**
 * Get Variant sites comparison among samples
 * @param {string} sample - Sample of interest
 * @param {number} position - Position in xAxis for tooltips display
 * @param {string} segment - Name of segment of interest
 * @param {WgsCovPlotDB} db
 * @returns {(string[])[]}
 */
function getSegmentVariantComparison(
    {
        sample,
        position,
        segment,
        db: {
            selectedSamples,
            depths,
            variants,
            segments_ref_id,
            segments_ref_seq,
        },
    }) {
    let variantArr = [];
    selectedSamples.forEach(sample => {
        let variantCall = find(variants[sample][segment],
            {"POS": position.toString()}, 0);
        if (!isNil(variantCall)) {
            variantArr.push(variantCall);
        } else {
            variantArr.push({
                "Sample": sample,
                "POS": position,
                "REF_ID": segments_ref_id[sample][segment],
                "Segment": segment,
                "Segment Length": depths[sample][segment].length,
                "REF_SEQ": segments_ref_seq[sample][segment][position - 1]
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
            selectedSamples.forEach(sample => {
                const segDepths = depths[sample][segment];
                let rowValue;
                if (segDepths.length === 0) {
                    rowValue = `No result reported for segment ${segment}`;
                } else {
                    if (position <= segDepths.length) {
                        rowValue = segDepths[position - 1].toLocaleString();
                    } else {
                        //Out of range when segment length < padding array
                        rowValue = `No sequence at this position. Reference sequence 
                            ${segments_ref_id[sample][segment]} is only 
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
                        row.push(variant.bold()); // Bold highlight selected sample
                    } else {
                        row.push(variant);
                    }
                } else {
                    row.push(""); // No information
                }
            });
        }
        rows.push(...[row]);
    });
    return rows;
}

/**
 * Get Coverage stat comparison among samples
 * @param {string[]} selectedSamples - List of sample names
 * @param {string} segment - Segment name
 * @param {Object} depths - Object of depths array
 * @param {number} position - Position in xAxis
 * @param {number} low - low coverage threshold
 * @returns {(string[])[]}
 */
function getSegmentCoverageStatComparison(
    {
        selectedSamples,
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
    selectedSamples.forEach((sample) => {
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
 * Get formatter for tooltips display
 * @param {Object} db - wgscovplot DB object
 * @returns {Array<Object>}
 */
function segmentTooltipFormatter({db}) {
    return function (params) {
        const {
            selectedSamples,
            depths,
            crossSampleComparisonInTooltips,
            variants,
            segments_ref_seq,
            segments_ref_id,
            low_coverage_threshold,
            showCovStatsInTooltips,
            segCoords
        } = db;
        let param = params[0];
        let output = "";
        // pos of xAxis in full scale
        let position = param.axisValue;
        let i = param.axisIndex;
        let positionRows = [];
        let coverageStatRows = [];
        let sample = selectedSamples[i];
        let segment = whichSegment({position, segCoords});
        let sequence = segments_ref_seq[sample][segment];
        let segmentLength = sequence.length;
        position = position - segCoords[segment].start + 1;
        let coverageDepth;
        let sampleSegDepths = depths[sample][segment];
        let refID = segments_ref_id[sample][segment];
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
        const isVariantBar = params.find(element => element.componentSubType === "bar");
        if (isVariantBar) { //tooltips at Variant sites
            if (crossSampleComparisonInTooltips) {
                positionRows = getSegmentVariantComparison(
                    {
                        sample,
                        position,
                        segment,
                        db
                    });
            } else {
                let foundObj = find(variants[sample][segment],
                    {"POS": position.toString()});
                if (!isNil(foundObj)) {
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
            if (showCovStatsInTooltips) {
                if (crossSampleComparisonInTooltips) {
                    coverageStatRows = getSegmentCoverageStatComparison(
                        {
                            selectedSamples,
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
            if (!isNil(db.primer_matches)) {
                let primerInfoRows = getPrimerInfo(sample, position , segment, db)
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
    };
}

/**
 * Get Variant sites comparison among samples
 * @param {string} sample - Sample of interest
 * @param {number} position - Position in xAxis for tooltips display
 * @param {string} segment - Name of segment of interest
 * @param {WgsCovPlotDB} db
 * @returns {(string[])[]}
 */
function getPrimerInfo (sample, position, segment, db){
    let primerInfo = [];
    if (has(db.primer_matches, [sample, segment])) {
        primerInfo = db.primer_matches[sample][segment];
    }
    let primerInfoRows = [];
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
    return primerInfoRows;
}

/**
 * Define options for tooltips
 * @param {Object} db - wgscovplot DB object
 * @returns {Array<Object>}
 */
function getSegmentTooltips(db) {
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
            formatter: segmentTooltipFormatter({db}),
        },
    ];
}

export {
    getSegmentTooltips
};
