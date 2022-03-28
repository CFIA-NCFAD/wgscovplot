import {slice} from "lodash/array";
import {filter, sortBy} from "lodash/collection";
import {mean} from "lodash/math";

/**
 * Calculate mean coverage for a genome region
 * @param {Array<number>} depths - depths array
 * @param {number} start - start position
 * @param {number} end - end position
 * @returns {number} - Returns the mean
 */
function meanCoverage(depths, start, end) {
    let subArray = slice(depths, start - 1, end);
    return mean(subArray);
}

/**
 *  Calculate coverage depth of genome region according to threshold low
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} low - the threshold that wants to set
 * @returns {number} - Returns the genome coverage
 */
function genomeCoverage(depths, start, end, low) {
    let subArray = slice(depths, start - 1, end);
    let filetedArray = filter(subArray, function (x) {
        return x >= low;
    });
    return (filetedArray.length / (end - start + 1)) * 100;
}

/**
 * Calculate median of an array
 * @param {Array<number>} arr - The array to iterate over
 * @returns {number} - Returns the median
 */
function median(arr) {
    let sortedArr = sortBy(arr);
    let half = Math.floor(sortedArr.length / 2);
    if (sortedArr.length % 2) {
        return sortedArr[half];
    } else{
        return (sortedArr[half - 1] + sortedArr[half]) / 2.0;
    }
}

/**
 * Calculate median coverage for a genome region
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @returns {number} - Returns the median
 */
function medianCoverage(depths, start, end) {
    let subArray = slice(depths, start - 1, end);
    return median(subArray);
}

export {median, meanCoverage, genomeCoverage, medianCoverage};