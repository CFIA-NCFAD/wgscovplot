import {ntColor} from "../../../util";
import {getSegmentsIndex} from "./getFluSegmentsInfo";
import {find} from "lodash/collection";

/**
 * Define variant sites bar
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Object} variants - The object of variants data
 * @param {Array<Array<number>>} segmentsRange - An array of segment start, end
 * @param {Object} refSeq - The object of reference sequences of each segment
 * @param {boolean} variantSites- whether to show variant sites information (tooltips options)
 * @param {boolean} showMutation - whether to show mutation below variant sites
 * @param {boolean} hideOverlapMutation - whether to hide overlapping mutation under variants sites
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 *  variants: { 'SAMPLE_NAME':{
 *                'SEGMENT_NAME': {}
 *              }
 *          }
 */
function getFluVariantSeries(samples, segments,
                             depths, variants,
                             segmentsRange, refSeq,
                             variantSites, showMutation,
                             hideOverlapMutation) {
    let variantSeries = [];
    let pos;
    for (let i = 0; i < samples.length; i++) {
        let data = []
        for (let j = 0; j < segments.length; j++) {
            for (let [k, varMap] of variants[samples[i]][segments[j]].entries()){
                pos = parseInt(varMap.POS);
                data.push([pos + segmentsRange[j][0] - 1, depths[samples[i]][segments[j]][pos - 1]]);
            }
        }
        variantSeries.push({
            type: "bar",
            xAxisIndex: i,
            yAxisIndex: i,
            data: data,
            barWidth: 2,
            itemStyle: {
                color: function (params) {
                    let axisValue = params.data[0];
                    let segmentIndex = getSegmentsIndex(axisValue, segmentsRange);
                    let nt = refSeq[samples[i]][segments[segmentIndex]][axisValue- segmentsRange[segmentIndex][0]];
                    if (ntColor.hasOwnProperty(nt)) {
                        return ntColor[nt];
                    }
                    return "#333";
                }
            },
            label: {
                show: showMutation,
                position: "bottom",
                align: "left",
                verticalAlign: "middle",
                distance: 10,
                color: "inherit",
                rotate: -30,
                formatter: function (params){
                    let output = "";
                    let segmentIndex = getSegmentsIndex(params.data[0], segmentsRange);
                    let foundObj = find(variants[samples[i]][segments[segmentIndex]],
                        {"POS": params.data[0] - segmentsRange[segmentIndex][0] + 1});
                    if (foundObj !== undefined && foundObj !== null){
                        output += foundObj.REF_SEQ + foundObj.POS + foundObj.ALT_SEQ
                    }
                    return output;
                }
            },
            labelLayout:{
                hideOverlap: hideOverlapMutation
            },
            large: true,
            tooltip:{
                trigger: variantSites ? "axis" : "none"
            }
        });
    }
    return variantSeries;
}


export {getFluVariantSeries};