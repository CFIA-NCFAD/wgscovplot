var _ = require(lodash)
const {median} = require('./wgscovplot.js');
var arr = [1, 4 , 7 , 9];

test('Median of an array [1, 4 , 7 , 9]', () => {
  expect(median(arr)).toBe(5.5);
});