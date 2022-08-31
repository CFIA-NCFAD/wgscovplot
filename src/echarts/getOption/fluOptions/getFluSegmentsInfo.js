import {segmentsColor} from "../../../util";
import {max, sum} from "lodash/math";

/**
 * Get maximum length for a segment among samples
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @returns {Array<number>}
 *
 */
function getMaxSegmentsLength (samples, segments, depths){
    let maxSegmentsLength = [];
    for (let i = 0; i < segments.length; i++){
        let maxLength = depths[samples[0]][segments[i]].length;
        for (let j = 0; j < samples.length; j++){
            if (maxLength <= depths[samples[j]][segments[i]].length){
                maxLength = depths[samples[j]][segments[i]].length;
            }
        }
        maxSegmentsLength.push(maxLength)
    }
    return  maxSegmentsLength;
}

/**
 * Get maxXAxis value
 * @param {Array<number>} maxSegmentsLength - An array of maximum length
 * @returns {number}
 */
function getXAxisMax (maxSegmentsLength){
    let xAxisMax = 0;
    xAxisMax = sum(maxSegmentsLength);
    return xAxisMax;
}

/**
 * Get maxYAxis value
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @returns {number}
 *
 */
function getYAxisMax (samples, segments, depths){
    let yAxisMax = 0;
    for (let i = 0; i < segments.length; i++){
        for (let j = 0; j < samples.length; j++){
            if (yAxisMax <= max(depths[samples[j]][segments[i]])){
                yAxisMax = max(depths[samples[j]][segments[i]])
            }
        }
    }
    return parseInt(yAxisMax * 1.5); // set max value for yAxis
}

/**
 * Get segment index which a position belongs to
 * @param {number} xAxisValue - xAxis value
 * @param {Array<Array<number>>} segmentsRange - An array segment start, end
 * @returns {number}
 */
function getSegmentsIndex (xAxisValue, segmentsRange){
    let index = 0;
    for (let i = 0; i < segmentsRange.length; i ++){
        if (xAxisValue >= segmentsRange[i][0] && xAxisValue <= segmentsRange[i][1]){
            index = i;
        }
    }
    return index;
}

/**
 * Get [start, end] for each segment
 * @param {Array<number>} maxSegmentsLength - An array of maximum length
 * @returns {Array<Array<number>>}
 *
 */
function getSegmentsRange (maxSegmentsLength){
    let segmentsRange = [];
    for (let m = 0; m < maxSegmentsLength.length; m++){
        let segmentsStart;
        let segmentsEnd;
        if (m == 0){
            segmentsStart = 1;
            segmentsEnd = maxSegmentsLength[m];
        }else{
            segmentsStart = segmentsRange.at(-1)[1] + 1;
            segmentsEnd = segmentsStart + maxSegmentsLength[m] - 1;
        }
        segmentsRange.push([segmentsStart, segmentsEnd]);
    }
    return segmentsRange;
}

/**
 * Get flu gene feature
 * @param {Array<string>} segments - An array of segments names
 * @param {Array<Array<number>>} segmentsRange - An array of maximum length
 * @returns {Array<Object>}
 *
 */
function getFluGeneFeature(segments, segmentsRange){
    let geneFeature = [];
    for (let i = 0; i < segments.length; i++){
        geneFeature.push({
            name: 'Segment ' + segments[i],
            value: {
                 "idx": i,
                 "start": segmentsRange[i][0],
                 "end": segmentsRange[i][1],
                 "level": 0,
                 "strand": '',
                 "rotate": 0.0,
                 "type": "segment_feature"
            },
            itemStyle:{
                "color": segmentsColor[segments[i]]
            }
        });
    }
    return geneFeature;
}
export {getMaxSegmentsLength, getSegmentsRange, getFluGeneFeature, getYAxisMax, getXAxisMax, getSegmentsIndex}