const {median, genomeCoverage, medianCoverage, meanCoverage} = require('../../wgscovplot/tmpl/js/wgscovplot.prod.bundle.js');
var depths = [
    [1, 2, 45, 3, 2, 34, 54, 65, 7, 6, 34, 45, 56, 67, 78, 78],
    [1, 2, 45, 0, 0, 9, 15, 65, 7, 6, 20, 8, 4, 15, 100, 102],
    [12, 12, 425, 3, 2, 10, 12, 9, 7, 6, 1, 45, 45, 67, 87, 97]
];

test('Median of an array ' + '[' + depths[0] + ']', () => {
    expect(median(depths[0])).toBe(39.5);
});

var gridIndex0 = 0;
test('Mean Coverage of depths array [' + depths[gridIndex0] + '], with region:' + 0 + '-' + (depths[gridIndex0].length - 3), () => {
    expect(meanCoverage(depths[gridIndex0], 1, depths[gridIndex0].length - 3).toFixed(2)).toBe("27.23");
});

var gridIndex1 = 1;
test('Genome Coverage of depths array [' + depths[gridIndex1] + '], with region:' + 0 + '-' + (depths[gridIndex1].length - 1), () => {
    expect(genomeCoverage(depths[gridIndex1], 1, depths[gridIndex1].length - 1, 10).toFixed(2)).toBe("40.00");
});

var gridIndex2 = 2;
test('Median Coverage of depths array [' + depths[gridIndex2] + '], with region:' + 1 + '-' + (depths[gridIndex2].length - 3), () => {
    expect(medianCoverage(depths[gridIndex2], 2, depths[gridIndex2].length - 3).toFixed(2)).toBe("9.50");
});

