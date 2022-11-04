import {FLU_SEGMENT_COLOURS} from "../../util";
import {max} from "lodash/math";


/**
 * Get maximum length for a segment among samples
 * @param {WgsCovPlotDB} db
 * @returns {SegmentCoords}
 */
function getSegmentCoords(db) {
    const {
        selectedSegments,
        selectedSamples,
        depths,
    } = db;
    let out = {};
    let prev = {
        maxLength: 0,
        start: 0,
        end: 1,
    }
    for (let segment of selectedSegments) {
        let maxLength = 0;
        for (let sample of selectedSamples) {
            let length = depths[sample][segment].length;
            if (maxLength <= length) {
                maxLength = length;
            }
        }
        let segCoords = {
            maxLength,
            start: prev.end + 1,
            end: prev.end + maxLength + 1,
        };
        out[segment] = segCoords;
        prev = segCoords;
    }
    return out;
}


/**
 * Get max y-axis value
 * @param {WgsCovPlotDB} db
 * @returns {number}
 */
function getYAxisMax(db) {
    const {
        selectedSamples,
        selectedSegments,
        depths
    } = db;
    let yAxisMax = 0;
    for (let segment of selectedSegments) {
        for (let sample of selectedSamples) {
            let maxDepth = max(depths[sample][segment]);
            if (yAxisMax <= maxDepth) {
                yAxisMax = maxDepth;
            }
        }
    }
    return yAxisMax;
}

/**
 * Find which segment a position lies in
 * @param {number} position - x-axis position value
 * @param {SegmentCoords} segCoords - Map of segment to segment coordinates
 * @returns {string | null} Name of segment if found
 */
function whichSegment(
    {
        position,
        segCoords
    }) {
    for (let segment of Object.keys(segCoords)) {
        let coords = segCoords[segment];
        if (position >= coords.start && position <= coords.end) {
            return segment
        }
    }
    return null;
}

/**
 * Get influenza segment features options array for ECharts
 * @param {WgsCovPlotDB} db
 * @returns {Array<Object>}
 *
 */
function getFluGeneFeature(db) {
    const {
        selectedSegments,
        segCoords
    } = db;
    let geneFeature = [];
    let i = 0
    for (let segment of selectedSegments) {
        let {start, end} = segCoords[segment];
        geneFeature.push({
            name: `Segment ${segment}`,
            value: {
                idx: i,
                start,
                end,
                level: 0,
                strand: "",
                rotate: 0.0,
                type: "segment",
            },
            itemStyle: {
                "color": FLU_SEGMENT_COLOURS[segment]
            }
        });
        i += 1;
    }
    return geneFeature;
}

export {getSegmentCoords, getFluGeneFeature, getYAxisMax, whichSegment};