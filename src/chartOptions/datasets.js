import {constant, times} from "lodash/util";
import {toFloat32Array} from "../util";

/**
 * Get ECharts dataset options containing depths for all samples (and segments if virus is segmented)
 * @param {WgsCovPlotDB} db
 * @returns {Object[]}
 */
function getDatasets(db) {
    const {
        selectedSamples,
        selectedSegments,
        depths,
        positions,
        segCoords,
    } = db
    let datasets = [];
    if (db.segment_virus === true) { // segment virus
        for (let sample of selectedSamples) {
            let depthArray = [];
            for (let segment of selectedSegments) {
                db.depths[sample][segment] = toFloat32Array(depths[sample][segment])
                let ds = db.depths[sample][segment];
                let coords = segCoords[segment];
                if (ds.length < coords.maxLength) {
                    // padding value 1E-5
                    let padding = times(coords.maxLength - ds.length, constant(1E-5));
                    ds = [...ds, ...padding]
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
            db.depths[sample] = toFloat32Array(depths[sample]);
            datasets.push({
                dimensions: [
                    {name: "depth", type: "float"},
                    {name: "position", type: "int"},
                ],
                source: {
                    position: positions,
                    depth: db.depths[sample],
                },
            });
        }
    }
    return datasets;
}

export {getDatasets};