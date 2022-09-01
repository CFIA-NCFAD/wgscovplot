import {times, constant} from "lodash/util";
import {concat} from "lodash/array";
import {getMaxSegmentsLength} from "./getFluSegmentsInfo";

/**
 * Get Flu Datasets
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @param {Array<number>} position - Array position for xAxis
 * @returns {Array<Object>}
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getFluDatasets(samples, segments, depths, position) {

    let datasets = [];
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
                position: position,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

export {getFluDatasets};