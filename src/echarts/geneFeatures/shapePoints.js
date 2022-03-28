/**
 * Define the points of gene/amplicon features shape
 * @param {number} x - x-axis coordinate
 * @param {number} y - y-axis coordinate
 * @param {number} width - width of shape
 * @param {number} height - height of shape
 * @param {number} strand - strand of feature
 * @param {string} feature - gene feature or amplicon feature
 * @returns {Array<Array<number>>} - Coordinate of 5 points for gene feature or 4 points for amplicon feature
 */
function shapePoints(x, y, width, height, strand, feature) {
    let taperOffset;
    if (width < 10) { // Element width is too small
        taperOffset = width/2;
    }else{
        taperOffset = 5;
    }
    if (feature === 'gene_feature') {
        if (strand === 1) {
            return [
                [x, y],
                [x + width - taperOffset, y],
                [x + width, y - height / 2],
                [x + width - taperOffset, y - height],
                [x, y - height],
            ];
        } else {
            return [
                [x, y - height / 2],
                [x + taperOffset, y],
                [x + width, y],
                [x + width, y - height],
                [x + taperOffset, y - height],
            ];
        }
    } else if (feature === 'amplicon_feature') {
        return [
            [x, y],
            [x + width, y],
            [x + width, y - height],
            [x, y - height],
        ];
    }
    else {
        return null;
    }
}

export {shapePoints};