import {times,constant} from "lodash/util";
import {concat} from "lodash/array";
import {getMaxSegmentsLength} from "./getFluSegmentsInfo";

/**
 * Get dataset information
 * @param {Array<Array<number>>} depths - Array of depths
 * @returns {Array<Object>}
 */
function getFluDatasets(samples, segments, depths, position) {

    let datasets = [];
    const maxSegmentsLength = getMaxSegmentsLength(samples, segments, depths)
    if (maxSegmentsLength.length == 0){
        return  datasets;
    }
    for (let k = 0; k < samples.length; k++) {
        let depthArray = [];
        for (let m=0; m <segments.length; m++){
            let temp = depths[samples[k]][segments[m]];
            if (temp.length < maxSegmentsLength[m]){
                let temp1 = times(maxSegmentsLength[m] - temp.length, constant(1E-10));
                temp = concat(depths[samples[k]][segments[m]], temp1)
            }
            depthArray = concat(depthArray, temp)
        }
        //console.log(samples[k], position.length, depthArray.length, depthArray);
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
    console.log(datasets)
    return datasets;
}

export {getFluDatasets};