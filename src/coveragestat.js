// Select methods for smaller webpack bundles
const slice = require('lodash/slice');
const filter = require('lodash/filter');
const sortBy = require('lodash/sortBy');
const mean = require('lodash/mean');

/**
 * Calculate mean coverage for a genome region
 * @param {Array<number>} depths - depths array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @returns {number}
 */
function meanCoverage(depths, start, end, gridIndex) {
    var subArray = slice(depths[gridIndex], start - 1, end);
    return mean(subArray);
}

/**
 *  Calculate genome coverage depth according to threshold low
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @param {number} low - the threshold that wants to set
 * @returns {number}
 */
function genomeCoverage(depths, start, end, gridIndex, low) {
    var subArray = slice(depths[gridIndex], start - 1, end);
    var filetedArray = filter(subArray, function (x) {
        return x >= low;
    });
    return (filetedArray.length / (end - start + 1)) * 100;
}

/**
 * Calculate median of an array
 * @param {Array<number>} arr
 * @returns {number}
 */
function median(arr) {
    var sortedArr = sortBy(arr)
    var half = Math.floor(sortedArr.length / 2);
    if (sortedArr.length % 2) return sortedArr[half];
    else return (sortedArr[half - 1] + sortedArr[half]) / 2.0;
}

/**
 * Calculate median coverage for a genome region
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @returns {number}
 */

function medianCoverage(depths, start, end, gridIndex) {
    var subArray = slice(depths[gridIndex], start - 1, end);
    return median(subArray);
}
export {median, meanCoverage, genomeCoverage, medianCoverage};