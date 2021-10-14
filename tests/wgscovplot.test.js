function median(arr) {
    arr.sort(function (a, b) {
        return a - b;
    });
    var half = Math.floor(arr.length / 2);
    if (arr.length % 2) return arr[half];
    else return (arr[half - 1] + arr[half]) / 2.0;
}

function sum(a, b) {
  return a + b;
}

var arr = [1, 4 , 7 , 9];

test('Median of an array [1, 4 , 7 , 9]', () => {
  expect(median(arr)).toBe(5.5);
});

test('adds 1 + 2 to equal 3', () => {
  expect(sum(1, 2)).toBe(3);
});