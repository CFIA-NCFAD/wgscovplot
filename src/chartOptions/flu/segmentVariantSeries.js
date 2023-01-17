import {NT_COLOURS} from "../../util";
import {whichSegment} from "./segmentInfo";
import {find} from "lodash/collection";
import {isEmpty} from "lodash/lang";
import {get, isNil} from "lodash";

/**
 * @param {string[]} segments - List of segments
 * @param {SegmentCoords} segCoords - An array of segment start, end
 * @returns {Object} - Coverage Chart Option
 */
function segmentSeparatorLines(segments, segCoords) {

    let data = [];
    for (let i = 0; i < segments.length; i++) {
        let segment = segments[i];
        if (i === 0) {
            data.push({xAxis: segCoords[segment].end});
        } else if (i === segments.length - 1) {
            data.push({xAxis: segCoords[segment].start});
        } else {
            data.push({xAxis: segCoords[segment].start});
            data.push({xAxis: segCoords[segment].end});
        }
    }
    return {
        silent: true,
        symbol: ["none", "none"],
        label: {
            show: false,
        },
        lineStyle: {
            color: "#000",
            width: 1,
            type: "dashed",
            opacity: 0.2
        },
        data: data
    };
}

/**
 * Define variant sites bar
 * @param {WgsCovPlotDB} db
 * @returns {Object} - Coverage Chart Option
 */
function getSegmentVariantSeries(db) {
    const {
        selectedSamples,
        selectedSegments,
        segCoords,
        variants,
        depths,
        hideOverlappingVariantLabels = true,
        showVariantSiteTooltips = true
    } = db;
    let variantSeries = [];
    let pos;
    for (let i = 0; i < selectedSamples.length; i++) {
        let data = [];
        let sample = selectedSamples[i];
        for (let segment of selectedSegments) {
            let vars = variants[sample][segment];
            if (!isNil(vars)) {
                for (let [k, varMap] of vars.entries()) {
                    pos = parseInt(varMap.POS);
                    data.push(
                        [
                            pos + segCoords[segment].start - 1,
                            depths[sample][segment][pos - 1]
                        ]);
                }
            } else {
                data.push([]);
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
                    const [position] = params.data;
                    const segment = whichSegment({position, segCoords});
                    const seqPosition = position - segCoords[segment].start;
                    const nt = db.segments_ref_seq[sample][segment][seqPosition];
                    return get(NT_COLOURS, nt, "#333");
                }
            },
            label: {
                show: db.showVariantLabels,
                position: "bottom",
                align: "left",
                verticalAlign: "middle",
                distance: 10,
                color: "inherit",
                rotate: -30,
                formatter: function (params) {
                    const [position] = params.data;
                    let segment = whichSegment({position, segCoords});
                    let seqPosition = position - segCoords[segment].start + 1;
                    let variant = find(variants[sample][segment], {POS: seqPosition.toString()}, 0);
                    if (!isNil(variant)) {
                        const {REF_SEQ, POS, ALT_SEQ} = variant;
                        return `${REF_SEQ}${POS}${ALT_SEQ}`;
                    }
                    return "";
                }
            },
            labelLayout: {
                hideOverlap: hideOverlappingVariantLabels
            },
            markLine: segmentSeparatorLines(selectedSegments, segCoords),
            large: true,
            tooltip: {
                trigger: showVariantSiteTooltips ? "axis" : "none"
            }
        });
    }
    return variantSeries;
}


export {getSegmentVariantSeries};