import {getSegmentsRange} from "./getFluSegmentsInfo";
import {toTableHtml} from "../../../util";


function getFluTooltips(samples, segments, depths, segmentsRange) {

    let toolTips = [
        {
            trigger: "axis",
            enterable: true,
            triggerOn: "mousemove",
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            position: "cursor",
            axisPointer: {
                type: 'line'
            },
            formatter: function (params) {
                let output = ''
                let rows = [];
                let coverageDepth = params[0].data[0];
                let pos = params[0].data[1];
                let segment = '';
                for (let i=0; i < segmentsRange.length; i ++){
                    if (pos >= segmentsRange[i][0] && pos <= segmentsRange[i][1]){
                        pos = pos - segmentsRange[i][0] + 1;
                        segment = segments[i];
                    }
                }
                if (coverageDepth == 1E-20){
                    coverageDepth = 0
                } else if (coverageDepth == 1E-10){
                    return 'Out of Range'
                }
                let segmentLength = depths[samples[params[0].axisIndex]][segment].length;
                output += "<h5>" + "Sample: " + samples[params[0].axisIndex] + "</h5>";
                rows = [
                    ['Segment', segment],
                    ['Segment Length', segmentLength.toLocaleString()],
                    ['Position', pos.toLocaleString()],
                    ['Coverage Depth', coverageDepth.toLocaleString()]
                ]
                output += toTableHtml(["Position Info", ""], rows, "table small");
                return output;
            }
        }
    ];
    return toolTips;
}

export {getFluTooltips};