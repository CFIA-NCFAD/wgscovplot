import {median} from '../wgscovplot/tmpl/js/wgscovplot.js';

function sum(a, b) {
  return a + b;
}

var arr = [1, 4 , 7 , 9];

test('Median of an array [1, 4 , 7 , 9]', () => {
  expect(median(arr)).toBe(5.5);
})

test('adds 1 + 2 to equal 3', () => {
  expect(sum(1, 2)).toBe(3);
});