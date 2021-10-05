function meanCoverage(depths, start, end, gridIndex) {
  var sub_array = _.slice(depths[gridIndex], start - 1, end);
  return _.mean(sub_array);
}

function genomeCoverage(depths, start, end, gridIndex, low) {
  var sub_array = _.slice(depths[gridIndex], start - 1, end);
  var fileted_array = _.filter(sub_array, function (x) {
    return x >= low;
  });
  return (fileted_array.length / (end - start + 1)) * 100;
}

function median(arr) {
  arr.sort(function (a, b) {
    return a - b;
  });
  var half = Math.floor(arr.length / 2);
  if (arr.length % 2) return arr[half];
  else return (arr[half - 1] + arr[half]) / 2.0;
}

function medianCoverage(depths, start, end, gridIndex) {
  var sub_array = _.slice(depths[gridIndex], start - 1, end);
  return median(sub_array);
}
