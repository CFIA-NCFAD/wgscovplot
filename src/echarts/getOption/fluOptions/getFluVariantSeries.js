import {ntColor} from "../../../util";
import {getSegmentsIndex} from "./getFluSegmentsInfo";


function getFluVariantSeries(samples, segments, depths, variants, segmentsRange, refSeq,
                             isVariantSites, isShowMutation, isHideOverlapMutation) {
    let variantSeries = [];
    let pos;
    for (let i = 0; i < samples.length; i++) {
        let data = []
        for (let j = 0; j < segments.length; j++) {
            for (let [k, varMap] of variants[samples[i]][segments[j]].entries()){
                pos = parseInt(varMap.POS);
                data.push([pos + segmentsRange[j][0] - 1, depths[samples[i]][segments[j]][pos -1]]);
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
                    let axisValue = params.data[0];
                    let segmentIndex = getSegmentsIndex(axisValue, segmentsRange);
                    let nt = refSeq[samples[i]][segments[segmentIndex]][axisValue- segmentsRange[segmentIndex][0]];
                    if (ntColor.hasOwnProperty(nt)) {
                        return ntColor[nt];
                    }
                    return "#333";
                }
            },
            large: true,
            tooltip:{
                trigger: isVariantSites ? "axis" : "none"
            }
        });
    }
    return variantSeries;
}


export {getFluVariantSeries};