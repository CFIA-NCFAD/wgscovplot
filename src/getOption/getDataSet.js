import {getMaxSegmentsLength} from "./fluOption/getFluSegmentsInfo";
import {constant, times} from "lodash/util";
import {concat} from "lodash/array";

/**
 * Get dataset information
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param depths - Array of depths, (depths is pobject for segments virus)
 * @param {Array<number>} positions - an array of genome positions which represent in X Axis
 * @returns {Array<Object>}
 *
 * segment virus depths format
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getDataSet(samples, segments, depths, positions) {
    let datasets = [];
    if (samples.length > 0 && segments.length > 0) { // segment virus
        const maxSegmentsLength = getMaxSegmentsLength(samples, segments, depths);
        for (let i = 0; i < samples.length; i++) {
            let depthArray = [];
            for (let j = 0; j < segments.length; j++) {
                let temp = depths[samples[i]][segments[j]];
                if (temp.length < maxSegmentsLength[j]) {
                    let temp1 = times(maxSegmentsLength[j] - temp.length, constant(1E-10)); // padding value 1E-10
                    temp = concat(depths[samples[i]][segments[j]], temp1);
                }
                depthArray = concat(depthArray, temp);
            }
            datasets.push({
                dimensions: [
                    {name: "depth", type: "float"},
                    {name: "position", type: "int"},
                ],
                source: {
                    position: positions,
                    depth: depthArray,
                },
            });
        }
    } else { // non-segment virus
        for (let [j, depthArray] of depths.entries()) {
            datasets.push({
                dimensions: [
                    {name: "depth", type: "float"},
                    {name: "position", type: "int"},
                ],
                source: {
                    position: positions,
                    depth: depthArray,
                },
            });
        }
    }
    return datasets;
}

export {getDataSet};