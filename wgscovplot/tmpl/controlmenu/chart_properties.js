function renderEnv() {
  var render_env = document.getElementById("renderenv").value;
  var isChecked = document.getElementById("toggledarkmode").checked;
  var samples = [];
  _.forEach(window.samples, function (value, key) {
    if (key < default_num_chart) {
      samples.push(value);
    }
  });
  var mode = "white";
  if (isChecked) mode = "dark";
  else mode = "white";
  if (render_env === "canvas") {
    echarts.dispose(chart); // destroy chart instance and re-init chart
    $chart = document.getElementById("chart");
    chart = echarts.init($chart, mode, { renderer: "canvas" });
    chart.setOption((option = getOption()));
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateChartOptionsMenu();
  } else {
    echarts.dispose(chart);
    $chart = document.getElementById("chart");
    chart = echarts.init($chart, mode, { renderer: "svg" });
    chart.setOption((option = getOption()));
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateChartOptionsMenu();
  }
}

$("#toggledarkmode").change(function () {
  var render_env = "canvas";
  if (document.getElementById("renderenv"))
    render_env = document.getElementById("renderenv").value;
  var samples = [];
  _.forEach(window.samples, function (value, key) {
    if (key < default_num_chart) {
      samples.push(value);
    }
  });
  if ($(this).prop("checked")) {
    echarts.dispose(chart); // destroy chart instance and re-init chart
    chart = echarts.init($chart, "dark", { renderer: render_env });
    chart.setOption((option = getOption()));
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateChartOptionsMenu();
  } else {
    echarts.dispose(chart); // destroy chart instance and re-init chart
    chart = echarts.init($chart, "white", { renderer: render_env });
    chart.setOption((option = getOption()));
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateChartOptionsMenu();
  }
});

function updateChartHeight(val) {
  document.getElementById("chartheightoutput").value = val + "%";
  var grid_option = chart.getOption().grid;
  for (var i = 0; i < grid_option.length; i++) {
    if (i < grid_option.length - 1)
      // don't change height of DNA feature charts
      grid_option[i]["height"] = val + "%";
    if (i > 0) {
      grid_option[i]["top"] =
        _.toInteger(grid_option[i - 1]["top"].replace("%", "")) +
        _.toInteger(grid_option[i - 1]["height"].replace("%", "")) +
        4 +
        "%";
    }
  }
  chart.setOption({ grid: grid_option });
}

function updateChartLeft(val) {
  document.getElementById("chartleftoutput").value = val + "%";
  var grid_option = chart.getOption().grid;
  _.forEach(grid_option, function (element) {
    element.left = val + "%";
  });
  chart.setOption({ grid: grid_option });
}

function updateChartRight(val) {
  document.getElementById("chartrightoutput").value = val + "%";
  var grid_option = chart.getOption().grid;
  _.forEach(grid_option, function (element) {
    element.right = val + "%";
  });
  chart.setOption({ grid: grid_option });
}

function updateChartTop(val) {
  document.getElementById("charttopoutput").value = val + "%";
  var grid_option = chart.getOption().grid;
  for (var i = 0; i < grid_option.length; i++) {
    if (i == 0) {
      grid_option[i]["top"] = val + "%";
    } else {
      grid_option[i]["top"] =
        _.toInteger(grid_option[i - 1]["top"].replace("%", "")) +
        _.toInteger(grid_option[i - 1]["height"].replace("%", "")) +
        _.toInteger(val) +
        "%";
    }
  }
  chart.setOption({ grid: grid_option });
}
