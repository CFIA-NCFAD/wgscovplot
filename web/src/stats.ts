import {mean, filter, sortBy, slice} from "lodash";

/**
 * Calculate mean coverage for a genome region
 * @param {Array<number>} depths - depths array
 * @param {number} start - start position
 * @param {number} end - end position
 * @returns {number} - Returns the mean
 */
function meanCoverage(depths: number[], start: number, end: number): number {
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
function genomeCoverage(depths: number[], start: number, end: number, low: number): number {
    let subArr = slice(depths, start - 1, end);
    let filteredArr = filter(subArr, function (x) {
        return x >= low;
    });
    return (filteredArr.length / (end - start + 1)) * 100;
}

/**
 * Calculate median of an array
 * @param {Array<number>} arr - The array to iterate over
 * @returns {number} - Returns the median
 */
function median(arr: number[]): number {
    let sortedArr = sortBy(arr);
    let half = Math.floor(sortedArr.length / 2);
    if (sortedArr.length % 2) {
        return sortedArr[half];
    } else {
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
function medianCoverage(depths: number[], start: number, end: number): number {
    let subArray = slice(depths, start - 1, end);
    return median(subArray);
}

export {median, meanCoverage, genomeCoverage, medianCoverage};