/**
 * Get Primer data
 * @param {string} sample
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} primerMatches - Object of primer data
 * @param {SegmentCoords} segCoords - An array of segment start, end
 * @returns {Object} - Object of primer data
 *
 */
import {graphic} from "echarts/core";
import {isEmpty} from "lodash/lang";

export function getFluPrimerData(
    {
        sample,
        db: {
            segments,
            primerMatches,
            segCoords,
        },
    }) {
    let primerFeature = [];
    for (let segment of segments) {
        let matchesToSampleSegment = primerMatches[sample][segment];
        let coords = segCoords[segment];
        for (let match of matchesToSampleSegment) {
            let {
                start,
                end,
                name,
                query_aligned,
                target_aligned,
                matched_aligned,
                cigar,
                edit_distance,
            } = match;
            // TODO: push object instead of array? does echarts need value to be an array instead of object?
            primerFeature.push({
                value: [
                    start + coords.start,
                    end + coords.start,
                    10, // height for primer annotation in log scale
                    name,
                    query_aligned,
                    target_aligned,
                    matched_aligned,
                    cigar,
                    edit_distance,
                ],
                itemStyle: {color: "#b71ae3"}
            });
        }
    }
    return primerFeature;
}

function primerMatchRenderer() {
    return function ({coordSys}, api) {
        let [startX, startY] = api.coord([api.value(0), api.value(2)]);
        let [endX, endY] = api.coord([api.value(1), 1]);
        let rectShape = graphic.clipRectByRect(
            {
                x: startX,
                y: startY,
                width: endX - startX,
                height: (endY - startY) * 0.3
            },
            coordSys
        );
        return rectShape && {
            type: "rect",
            shape: rectShape,
            style: api.style(),
        };
    };
}

/**
 * Define primers series, it is used to annotate into coverage plot
 * @param {Object} db - wgscovplot DB object
 * @returns {Object} - Primer Series Option
 *
 */
function getFluPrimerSeries(db) {
    const {
        selectedSamples,
        primerMatches,
    } = db;
    let primerSeries = [];
    if (isEmpty(primerMatches)) {
        return primerSeries;
    }
    for (let [i, sample] of selectedSamples.entries()) {
        primerSeries.push({
            type: "custom",
            xAxisIndex: i,
            yAxisIndex: i,
            renderItem: primerMatchRenderer(),
            label: {
                show: true,
                position: "top",
                distance: 20,
                rotate: 45
            },
            encode: {
                x: [0, 1],
                y: 3,
            },
            silent: true,
            data: getFluPrimerData({sample, db}),
        });
    }
    return primerSeries;
}

export {getFluPrimerSeries};