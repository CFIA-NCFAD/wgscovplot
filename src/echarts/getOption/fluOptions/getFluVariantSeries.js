import {ntColor} from "../../../util";
import {getSegmentsIndex} from "./getFluSegmentsInfo";


function getFluVariantSeries(samples, segments, depths, variants, segmentsRange, refSeq) {
    let variantSeries = [];
    let pos;
    for (let i = 0; i < samples.length; i++) {
        let data = []
        for (let j = 0; j < segments.length; j++) {
            for (let [k, varMap] of variants[samples[i]][segments[j]].entries()){
                pos = parseInt(varMap.POS);
                data.push([pos + segmentsRange[j][0] - 1, depths[samples[i]][segments[j]][pos]]);
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
                    let segmentIndex = getSegmentsIndex(params.data[0], segmentsRange);
                    let nt = refSeq[samples[i]][segments[segmentIndex]][params.data[0] - segmentsRange[segmentIndex][0]];
                    if (ntColor.hasOwnProperty(nt)) {
                        return ntColor[nt];
                    }
                    return "#333";
                }
            },
            large: true,
        });
    }
    return variantSeries;
}


export {getFluVariantSeries};