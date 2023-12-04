import {graphic} from "echarts/core";
import {isEmpty, has, isNil} from "lodash";
import {WgsCovPlotDB} from "../../db";

export function getSegmentPrimerData(sample: string, db: WgsCovPlotDB) {
    let primerFeatures = [];
    let segments = Object.keys(db.segCoords); // trigger only segCoords is updated from get Dataset
    for (let segment of segments) {
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

export function primerMatchRenderer() {
    // @ts-ignore
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

export function getSegmentPrimerSeries(db: WgsCovPlotDB) {
    console.log("Trigger getSegmentPrimerSeries")
    let primerSeries: any[] = [];
    if (isEmpty(db.primer_matches)) {
        return primerSeries;
    }
    for (let [i, sample] of db.chartOptions.selectedSamples.entries()) {
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
            data: getSegmentPrimerData(sample, db),
        });
    }
    //console.timeEnd("Trigger getSegmentPrimerSeries")
    return primerSeries;
}

export function getPrimerInfo (sample: string, position: number, segment: string, db: WgsCovPlotDB){
    let primerInfo: any[] = [];
    if (has(db.primer_matches, [sample, segment])) {
        primerInfo = db.primer_matches[sample][segment];
    }
    let primerInfoRows: any[] = [];
    for (let i = 0; i < primerInfo.length; i++) {
        let {
            start,
            end,
            name,
            query_aligned,
            matched_aligned,
            target_aligned,
            edit_distance,
            cigar,
            other_locations,
        } = primerInfo[i];
        if (position >= start && position <= end) {
            primerInfoRows = [
                ["Primer Name", name],
                ["Primer Sequence", query_aligned],
                ["Match Aligned", matched_aligned],
                ["Ref Sequence", target_aligned],
                ["Cigar", cigar],
                ["Start", start + 1],
                ["End", end + 1],
                ["Edit Distance", edit_distance],
                ["Other Locations", other_locations]
            ];
        }
    }
    return primerInfoRows;
}