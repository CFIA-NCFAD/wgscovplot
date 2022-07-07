import {toTableHtml} from "../../../util";
import {getMaxSegmentsLength, getSegmentsIndex, getSegmentsRange} from "./getFluSegmentsInfo";
import {find} from "lodash/collection";
import {genomeCoverage, meanCoverage, medianCoverage} from "../../../coverageStat";

/**
 * Get Variant sites comparison among samples
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} variants - Object of variant
 * @param {Object} depths - Object of depths array
 * @param {Object} refSeq - The object of reference sequences of each segment
 * @param {Object} refID - The object of reference ID of each segment
 * @param {number} position - Position in xAxis
 * @param {number} segmentIndex - Segment Index
 * @param {string} sample - Selected sample
 * @returns {Array<Array<string>>}
 */
function getFluVariantComparison(samples, segments,
                                 depths, variants,
                                 refSeq, refID,
                                 position, segmentIndex, sample) {
    let rows = [];
    let variantArr = [];
    for (let i = 0; i < samples.length; i++){
        let foundObj = find(variants[samples[i]][segments[segmentIndex]],
        {"POS": position});
        if (foundObj !== undefined && foundObj !== null){
            variantArr.push(foundObj);
        } else {
            variantArr.push({"Sample": samples[i],
                             "POS": position,
                             "REF_ID": refID[samples[i]][segments[segmentIndex]],
                             "Segment": segments[segmentIndex],
                             "Segment Length": depths[samples[i]][segments[segmentIndex]].length,
                             "REF_SEQ": refSeq[samples[i]][segments[segmentIndex]][position-1]
                            });
        }
    }
    let unionKeys = [...new Set(variantArr.reduce((r, e) => [...r, ...Object.keys(e)], []))];
    unionKeys.push("Coverage Depth"); // Add Coverage Depth row
    unionKeys.forEach(key => {
        let row = [];
        row.push(key);
        if (key === "Coverage Depth"){
            for (let j = 0; j < samples.length; j++){
                let depth  = depths[samples[j]][segments[segmentIndex]];
                if (depth.length == 0){
                    row.push(`No result reported for segment ${segments[segmentIndex]}`);
                }else {
                    if (position <= depth.length){
                        row.push(depths[samples[j]][segments[segmentIndex]][position-1].toLocaleString());
                    }else{ //Out of range when segment length < padding array
                        row.push(`No sequence at this position. Reference sequence 
                            ${refID[samples[j]][segments[segmentIndex]]} is only 
                            ${depths[samples[j]][segments[segmentIndex]].length} bp`);
                    }
                }
            }
        }else {
            variantArr.forEach(element => {
                if (element[key] !== undefined && element[key] !== null){
                    if (key === "Sample" && element[key] === sample) {// Bold highlight selected sample
                        row.push(element[key].bold());
                    }else{
                        row.push(element[key]);
                    }
                }else {
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
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {number} position - Position in xAxis
 * @param {number} segmentIndex - Segment Index
 * @returns {Array<Array<string>>}
 */
function getFluCoverageStatComparison (samples, segments, depths, position, segmentIndex){
    let rows = [];
    let tableHeader = ["Sample", "Depth at position "+ position.toLocaleString(), "Segment",
        "Range", "Segment Length", "Mean Coverage (X)", "Median Coverage (X)", "Genome Coverage (>=10X) (%)"];
    rows.push(...[tableHeader]);
    for (let [i, sample] of samples.entries()){
        let depthArr = depths[samples[i]][segments[segmentIndex]]
        let meanCov = meanCoverage(depthArr, 1, depthArr.length).toFixed(2);
        let medianCov = medianCoverage(depthArr, 1, depthArr.length).toFixed(2);
        let genomeCov = genomeCoverage(depthArr, 1, depthArr.length, 10).toFixed(2);
        let coverageDepth
        if (position <= depthArr.length){
            coverageDepth = depthArr[position-1].toLocaleString();
        } else if (depthArr.length == 0){
            coverageDepth = `No result reported for segment ${segments[segmentIndex]}`
        }
        else {
            coverageDepth = "No sequence at this position"
        }
        let row = [sample, coverageDepth, segments[segmentIndex] ,1 + " - " + depthArr.length.toLocaleString(),
            depthArr.length.toLocaleString(), meanCov, medianCov, genomeCov];
        rows.push(...[row]);
    }
    return rows;
}

/**
 * Get tooltips options
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Object} variants - The object of variants data
 * @param {Object} refSeq - The object of reference sequences of each segment
 * @param {Object} refID - The object of reference ID of each segment
 * @param {string} triggerOnType - mousemove or click (tooltips options)
 * @param {boolean} coverageStatView - whether to show coverage statistics (tooltips options)
 * @param {boolean} infoComparison - whether to compare variants/coverage stat across samples (tooltips options)
 * @returns {Array<Object>} - tooltips Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 *  variants: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': {}
 *              }
 *          }
 *  refSeq: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': ref_seq
 *              }
 *          }
 *   refID: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': ref_id
 *              }
 *          }
 */
function getFluTooltips(samples, segments, depths, variants, refSeq, refID,
                        triggerOnType, infoComparison, coverageStatView) {

    let maxSegmentsLength = getMaxSegmentsLength(samples, segments, depths);
    let segmentsRange = getSegmentsRange(maxSegmentsLength);
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
                let param = params[0]
                let output = ''
                let position = param.axisValue; // pos of xAxis in full scale
                let i = param.axisIndex;
                let positionRows = [];
                let coverageStatRows = [];
                let sample = samples[i]
                let segmentIndex = getSegmentsIndex(position, segmentsRange); // find segment index which pos belongs to
                let segmentName = segments[segmentIndex];
                let segmentLength = refSeq[sample][segments[segmentIndex]].length; // get segment length
                position = position - segmentsRange[segmentIndex][0] + 1; // convert to pos in segment
                let coverageDepth;
                if (segmentLength == 0){
                    coverageDepth = `No result reported for segment ${segments[segmentIndex]}`
                }else{
                    if (position <= segmentLength) {
                        coverageDepth = depths[sample][segments[segmentIndex]][position-1].toLocaleString(); // get coverage depth for pos
                    }else{
                        coverageDepth = `No sequence at this position. 
                        Reference sequence ${refID[samples[i]][segments[segmentIndex]]} 
                        is only ${depths[samples[i]][segments[segmentIndex]].length} bp`
                    }
                }
                const variantBar = params.find(element => { // check if current position is variant sites
                    if (element.componentSubType === "bar") {
                        return true;
                    }
                    return false;
                });
                if (variantBar){ //tooltips at Variant sites
                    if (infoComparison){
                        positionRows = getFluVariantComparison(samples, segments, depths, variants,
                            refSeq, refID, position, segmentIndex , sample);
                    }else {
                        let foundObj = find(variants[samples[param.axisIndex]][segments[segmentIndex]],
                            {"POS": position});
                        if (foundObj !== undefined && foundObj !== null){
                            for (const [key, value] of Object.entries(foundObj)) {
                                if (key !== 'Sample') { // do not write row for sample name
                                    positionRows.push(
                                        ...[[key, value]]
                                    );
                                }
                            }
                        }
                        positionRows.push(["Coverage Depth", coverageDepth])
                    }
                } else { // tooltips for Non-Variant Sites
                    if (position > segmentLength){ //Out of range when segment length < padding array
                        positionRows = [
                            ['Position', position],
                            [coverageDepth, '']
                        ]
                        output += "<h5>" + "Sample: " + samples[params[0].axisIndex] + "</h5>";
                        output += toTableHtml(["Position Info", ""], positionRows, "table small");
                        return output
                    } // Pos within segment length
                    positionRows = [
                        ['Segment', segmentName],
                        ['Segment Length', segmentLength.toLocaleString()],
                        ['REF_ID', refID[samples[param.axisIndex]][segments[segmentIndex]]],
                        ['POS', position.toLocaleString()],
                        ['REF_SEQ', refSeq[samples[param.axisIndex]][segments[segmentIndex]][position-1]],
                        ['ALT_SEQ',''],
                        ['ALT_FREQ',''],
                        ['Coverage Depth', coverageDepth]
                    ]
                }
                if (positionRows.length){ // write rows to table
                    output += "<h5>" + "Sample: " + samples[params[0].axisIndex] + "</h5>";
                    output += toTableHtml(["Position Info", ""], positionRows, "table small");
                    if (coverageStatView){
                        if (infoComparison){
                            coverageStatRows = getFluCoverageStatComparison (samples, segments, depths, position, segmentIndex);
                        }
                        else{
                            let meanCov = meanCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength).toFixed(2);
                            let medianCov = medianCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength).toFixed(2);
                            let genomeCov = genomeCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength, 10).toFixed(2);
                            coverageStatRows = [
                                ["Mean Coverage", meanCov + "X"],
                                ["Median Coverage", medianCov + "X"],
                                ["Genome Coverage ( >= 10X)", genomeCov + "%"],
                            ];
                        }
                        output += toTableHtml(["Coverage View Stats", ""], coverageStatRows, "table small");
                    }
                    return output;
                }
            }
        }
    ];
    return toolTips;
}

export {getFluTooltips};