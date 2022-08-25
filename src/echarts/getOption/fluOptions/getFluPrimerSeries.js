import {graphic} from "echarts/core";
import {isEmpty} from "lodash/lang";

function getFluPrimerData(sample, segments, primerData, segmentsRanges){
    let primerFeature = [];
    for (let i = 0; i < segments.length; i++){
        let temp = primerData[sample][segments[i]]
        for(let j = 0; j < temp.length; j++){
            primerFeature.push({
                value: [temp[j]['start'] + segmentsRanges[i][0],
                        temp[j]['end'] + segmentsRanges[i][0],
                        10,
                        temp[j]['name'],
                        temp[j]['query_aligned'],
                        temp[j]['target_aligned'],
                        temp[j]['matched_aligned'],
                        temp[j]['cigar'],
                        temp[j]['edit_distance']],
                itemStyle: {"color": "#b71ae3"}
            });
        }
    }
    return primerFeature;
}

function getFluPrimerSeries(samples, segments, primerData, segmentsRanges) {
    let primerSeries = [];
    if (isEmpty(primerData)){
        return primerSeries;
    }
    for (let [i, sample] of samples.entries()) {
        let data = getFluPrimerData(sample, segments, primerData, segmentsRanges);
        primerSeries.push({
            type: "custom",
            xAxisIndex: i,
            yAxisIndex: i,
            renderItem: function (params, api) {
                let start = api.coord([api.value(0), api.value(2)]);
                let end = api.coord([api.value(1), 1]);
                let rectShape = graphic.clipRectByRect(
                    {
                        x: start[0],
                        y: start[1],
                        width: end[0] - start[0],
                        height: (end[1] - start[1]) * 0.3
                    },
                    {
                        x: params.coordSys.x,
                        y: params.coordSys.y,
                        width: params.coordSys.width,
                        height: params.coordSys.height
                    }
                );
                return rectShape && {
                    type: "rect",
                    shape: rectShape,
                    style: api.style(),
                };
            },
            label: {
                show: true,
                position: "top",
                distance: 20,
                rotate:45
            },
            encode: {
                x: [0, 1],
                y: 3,
            },
            silent: true,
            data: data,
        });
    }
    return primerSeries;
}

export {getFluPrimerSeries}