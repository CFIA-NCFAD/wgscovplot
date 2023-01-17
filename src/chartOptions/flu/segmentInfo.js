import {FLU_SEGMENT_COLOURS} from "../../util";

/**
 * Get maximum length for a segment among samples
 * @param {WgsCovPlotDB} db
 * @returns {Object}
 *
 */
function getMaxSegmentLength(db) {
    const {
        selectedSegments,
        selectedSamples
    } = db;
    let out = {}
    for (let segment of selectedSegments) {
        let maxLength = 0;
        for (let sample of selectedSamples) {
            let length = db.mosdepth_info[sample][segment].ref_seq_length;
            if (maxLength <= length) {
                maxLength = length;
            }
        }
        out[segment] = maxLength;
    }
    return out
}

/**
 * Get maximum length for a segment among samples
 * @param {WgsCovPlotDB} db
 * @param {maxSegmentLength} maxSegmentLength
 * @returns {Object}
 */
function getSegmentCoords(db, maxSegmentLength) {
    const {
        selectedSegments,
    } = db;
    let out = {};
    let prev = {
        maxLength: 0,
        start: 0,
        end: 0,
    }
    for (let segment of selectedSegments) {
        let maxLength = maxSegmentLength[segment];
        let segCoords = {
            maxLength,
            start: prev.end + 1,
            end: prev.end + maxLength,
        };
        out[segment] = segCoords;
        prev = segCoords;
    }
    return out;
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
 * @returns {Object[]}
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
                strand: 1,
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

export {getMaxSegmentLength, getSegmentCoords, getFluGeneFeature, whichSegment};