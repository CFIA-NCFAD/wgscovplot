import {getVariantComparison, toTableHtml} from "../../../util";
import {getSegmentsIndex} from "./getFluSegmentsInfo";
import {find} from "lodash/collection";
import {genomeCoverage, meanCoverage, medianCoverage} from "../../../coverageStat";

function getFluVariantComparison(samples, segments, variants, depths, segmentIndex, position, sample) {
    let rows = [];
    let variantArr = [];
    for (let i = 0; i < samples.length; i++){
        let foundObj = find(variants[samples[i]][segments[segmentIndex]],
        {"POS": position});
        if (foundObj !== undefined && foundObj !== null){
            variantArr.push(foundObj);
        } else {
            variantArr.push({"sample": samples[i], "POS":position});
        }
    }
    var unionKeys = [...new Set(variantArr.reduce((r, e) => [...r, ...Object.keys(e)], []))];
    unionKeys.push("Coverage Depth"); // Add Coverage Depth row
    unionKeys.forEach(key => {
        let row = [];
        row.push(key);
        if (key === "Coverage Depth"){
            for (let j = 0; j < samples.length; j++){
                row.push(depths[samples[j]][segments[segmentIndex]][position-1].toLocaleString());
            }
        }else {
            variantArr.forEach(element => {
                if (element[key] !== undefined && element[key] !== null){
                    if (key === 'sample' && element[key] === sample) {// Bold highlight selected sample
                        row.push(element[key].bold());
                    }else{
                        row.push(element[key]);
                    }
                }else {
                    row.push("");
                }
            });
        }
        rows.push(...[row]);
    });
    return rows;
}

function getFluTooltips(samples, segments, depths, variants, refSeq, segmentsRange,
                        triggerOnType, isInfoComparison, isCovergateStatView) {

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
                let sample = samples[i]
                let segmentIndex = getSegmentsIndex(position, segmentsRange); // find segment index which pos belongs to
                let segmentName = segments[segmentIndex];
                position = position - segmentsRange[segmentIndex][0] + 1; // covert to pos in segment
                let coverageDepth = depths[sample][segments[segmentIndex]][position-1]; // get coverage depth for pos
                let segmentLength = refSeq[sample][segments[segmentIndex]].length; // get segment length
                if (coverageDepth == 1E-10){
                    return 'Out of Range'
                }
                let positionRows = [];
                let coverageStatRows = [];
                const isVariantBar = params.find(element => {
                    if (element.componentSubType === "bar") {
                        return true;
                    }
                    return false;
                });
                if (isVariantBar){
                    if (isInfoComparison){
                        positionRows = getFluVariantComparison(samples, segments, variants, depths, segmentIndex, position, sample);
                    }else {
                        positionRows = [
                            ["Position", position.toLocaleString()],
                            ["Coverage Depth", coverageDepth.toLocaleString()],
                        ];
                        let foundObj = find(variants[samples[param.axisIndex]][segments[segmentIndex]],
                            {"POS": position});
                        if (foundObj !== undefined && foundObj !== null){
                            for (const [key, value] of Object.entries(foundObj)) {
                                if (key !== 'POS' && key !== 'sample') {
                                    positionRows.push(
                                        ...[[key, value]]
                                    );
                                }
                            }
                        }
                        positionRows.splice(3,0, ['Segment Length', segmentLength.toLocaleString()])
                    }
                } else {
                    positionRows = [
                        ['Position', position.toLocaleString()],
                        ['Coverage Depth', coverageDepth.toLocaleString()],
                        ['Segment', segmentName],
                        ['Segment Length', segmentLength.toLocaleString()],
                        ['REF', refSeq[samples[param.axisIndex]][segments[segmentIndex]][position-1]],
                    ]
                }
                if (positionRows.length){
                    output += "<h5>" + "Sample: " + samples[params[0].axisIndex] + "</h5>";
                    output += toTableHtml(["Position Info", ""], positionRows, "table small");
                    if (isCovergateStatView){
                        if (!isInfoComparison){
                            coverageStatRows = "";
                        }
                        else{
                            let meanCov = meanCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength).toFixed(2);
                            let medianCov = medianCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength).toFixed(2);
                            let genomeCov = genomeCoverage(depths[sample][segments[segmentIndex]], 1, segmentLength, 10).toFixed(2);
                            coverageStatRows = [
                                ["Mean Coverage", meanCov + "X"],
                                ["Median Coverage", medianCov + "X"],
                                ["Genome Coverage ( >= 10x)", genomeCov + "%"],
                            ];
                            output += toTableHtml(["Coverage View Stats", ""], coverageStatRows, "table small");
                        }
                    }
                    return output;
                }
            }
        }
    ];
    return toolTips;
}

export {getFluTooltips};