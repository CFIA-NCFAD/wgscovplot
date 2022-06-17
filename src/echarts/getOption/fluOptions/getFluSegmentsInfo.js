import {segmentsColor} from "../../../util";
import {max} from "lodash/math";

/**
 * Get maximum length for each segment
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of segments name
 * @param {Object} depths - Object of depths array
 * @returns {Array<number>}
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getMaxSegmentsLength (samples, segments, depths){
    let maxSegmentsLength = [];
    for (let i = 0; i < segments.length; i++){
        let maxLength = depths[samples[0]][segments[i]].length;
        for (let j = 0; j < samples.length; j++){
            // Find max depths.length for each segment
            if (maxLength <= depths[samples[j]][segments[i]].length){
                maxLength = depths[samples[j]][segments[i]].length;
            }
        }
        maxSegmentsLength.push(maxLength)
    }
    return  maxSegmentsLength;
}

/**
 * Get max depth among samples/segments
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of segments name
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
    return yAxisMax * 1.5; // set max value for yAxis
}

/**
 * Get range [start, end] for each segment
 * @param {Array<number>} maxSegmentsLength - An array of maximum length
 * @returns {Array<Array<number>>}
 *
 */
function getSegmentsRange (maxSegmentsLength){
    let segmentsRange = [];
    for (let m = 0; m < maxSegmentsLength.length; m++){
        let segmentsStart;
        let segmentsEnd;
        if (m==0){
            segmentsStart = 1;
            segmentsEnd = maxSegmentsLength[m];
        }else{
            segmentsStart = 1
            for (let n = 0; n < m ; n++){
                segmentsStart =  segmentsStart + maxSegmentsLength[n];
            }
            segmentsEnd = segmentsStart + maxSegmentsLength[m] - 1;
        }
        segmentsRange.push([segmentsStart, segmentsEnd]);
    }
    return segmentsRange;
}

/**
 * Get flu gene feature
 * @param {Array<string>} segments - An array of segments name
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
                 "strand": 1,
                 "rotate": 0.0,
                 "type": "gene_feature"
            },
            itemStyle:{
                "color": segmentsColor[segments[i]]
            }
        });
    }
    return geneFeature;
}
export {getMaxSegmentsLength, getSegmentsRange, getFluGeneFeature, getYAxisMax}