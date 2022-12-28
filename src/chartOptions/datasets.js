import {constant, times} from "lodash/util";
import {isNil} from "lodash";
import {toFloat32Array} from "../util";

/**
 * Get ECharts dataset options containing depths for all samples (and segments if virus is segmented)
 * @param {Array<string>} selectedSamples - An array of samples names
 * @param {SampleSegmentDepths | SampleDepths} selectedDepths
 * @param {Array<number>} positions - an array of genome positions which represent in X Axis
 * @param {Array<string> | null} segments - An array of segments names
 * @param {SegmentCoords | null} segCoords
 * @returns {Array<Object>}
 */
function getDatasets(db) {
    const {
        selectedSamples,
        depths,
        positions,
        segments= null,
        segCoords= null,
    } = db
    let datasets = [];
    if (!isNil(segments) && !isNil(segCoords)) { // segment virus
        for (let sample of selectedSamples) {
            let depthArray = [];
            for (let segment of segments) {
                let ds = depths[sample][segment];
                let coords = segCoords[segment];
                if (ds.length < coords.maxLength) {
                    // padding value 1E-10
                    let padding = times(coords - ds.length, constant(1E-10));
                    ds = [...ds, ...padding];
                }
                depthArray = [...depthArray, ...ds];
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
        for (let sample of selectedSamples) {
            datasets.push({
                dimensions: [
                    {name: "depth", type: "float"},
                    {name: "position", type: "int"},
                ],
                source: {
                    position: positions,
                    depth: toFloat32Array(depths[sample]),
                },
            });
        }
    }
    return datasets;
}

export {getDatasets};