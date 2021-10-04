function setScale() {
  var scale_type = document.getElementById("scale").value;
  var yaxis_option = chart.getOption().yAxis;
  if (scale_type == "value") {
    _.forEach(yaxis_option, function (element) {
      if (element.gridIndex < yaxis_option.length - 1) {
        element.type = scale_type;
        element.min = 0;
        element.max = 40000;
      }
    });
    chart.setOption({ yAxis: yaxis_option });
    document.getElementById("ymax").value = 40000;
  } else {
    _.forEach(yaxis_option, function (element) {
      if (element.gridIndex < yaxis_option.length - 1) {
        element.type = scale_type;
        element.min = 1;
        element.max = 100000;
      }
    });
    chart.setOption({ yAxis: yaxis_option });
    document.getElementById("ymax").value = 100000;
  }
}
function setYMax() {
  var ymax = document.getElementById("ymax").value;
  var yaxis_option = chart.getOption().yAxis;
  _.forEach(yaxis_option, function (element) {
    if (element.gridIndex < yaxis_option.length - 1) {
      element.max = ymax;
    }
  });
  chart.setOption({ yAxis: yaxis_option });
}
