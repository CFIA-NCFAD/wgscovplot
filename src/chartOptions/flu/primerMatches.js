
import {graphic} from "echarts/core";
import {isNil, has} from "lodash";

/**
 * Get Primer data
 * @param {string} sample
 * @param {WgsCovPlotDB} db
 * @returns {Object} - Object of primer data
 *
 */
function getSegmentPrimerData({sample, db}) {
    let primerFeatures = [];
    for (let segment of db.selectedSegments) {
        let matchesData;
        if (has(db.primer_matches, [sample, segment])) {
            matchesData = db.primer_matches[sample][segment];
        }
        let coords = db.segCoords[segment]
        if (!isNil(matchesData)) {
            for (let match of matchesData) {
                let {
                    start,
                    end,
                    name,
                    query_aligned,
                    target_aligned,
                    matched_aligned,
                    cigar,
                    edit_distance,
                    other_locations,
                } = match;
                // TODO: push object instead of array? does echarts need value to be an array instead of object?
                primerFeatures.push({
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
                        other_locations
                    ],
                    itemStyle: {color: "#b71ae3"}
                });
            }
        }
    }
    return primerFeatures;
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
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @returns {Object} - Primer Series Option
 *
 */
function getSegmentPrimerSeries(db) {
    const {
        selectedSamples,
    } = db;
    let primerSeries = [];
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
            data: getSegmentPrimerData({sample, db}),
        });
    }
    return primerSeries;
}

export {getSegmentPrimerSeries};