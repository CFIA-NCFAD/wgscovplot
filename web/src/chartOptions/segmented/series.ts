import {whichSegment} from "./getFluSegmentsInfo";
import {isEmpty, find} from "lodash";
import {remove, map, get, isNil} from "lodash";

/**
 * @typedef {Object<string, Object<string, Object<string, string>[]>>} SampleSegmentVariants
 */

/**
 * @param {SegmentCoords} segmentsInterval - An array of segment start, end
 */
function segmentSeparatorLines(segmentsInterval) {
    const starts = remove(map(segmentsInterval, "start"), (x) => x <= 1)
    const data = map(starts, function (x) {
        return {xAxis: x};
    });
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
function series(db) {
    const {
        selectedSamples,
        selectedSegments,
        segCoords,
        variants,
        depths,
    } = db;
    let variantSeries = [];
    let pos;
    for (let i = 0; i < selectedSamples.length; i++) {
        let data = [];
        let sample = selectedSamples[i];
        for (let segment of selectedSegments) {
            let vars = variants[sample][segment];
            if (!isNil(vars) || !isEmpty(vars)) {
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
                    const nt = db.ref_seqs[sample][segment][seqPosition];
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
                    let variant = find(variants[sample][segment], {POS: seqPosition});
                    if (!isNil(variant)) {
                        const {REF_SEQ, POS, ALT_SEQ} = variant;
                        return  `${REF_SEQ}${POS}${ALT_SEQ}`;
                    }
                    return "";
                }
            },
            labelLayout: {
                hideOverlap: hideOverlapMutation
            },
            markLine: segmentSeparatorLines(segCoords),
            large: true,
            tooltip: {
                trigger: variantSites ? "axis" : "none"
            }
        });
    }
    return variantSeries;
}


export {series};