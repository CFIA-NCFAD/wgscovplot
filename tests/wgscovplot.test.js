const lodash = require('lodash')
module.exports = {
    globals: {
        '_': lodash
    }
};
import {median} from "../src/coveragestat.js";
var arr = [1, 4 , 7 , 9];

test('Median of an array [1, 4 , 7 , 9]', () => {
  expect(median(arr)).toBe(5.5);
});