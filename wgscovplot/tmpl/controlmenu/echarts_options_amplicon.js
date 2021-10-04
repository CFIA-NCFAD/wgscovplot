function renderGeneFeatures(params, api) {
  var categoryIndex = api.value(0);
  var points;
  var shape;
  var style;
  var textConfig;
  var textContent;
  var height = gene_features_properties["rec_items_height"];
  var start = api.coord([api.value(1), categoryIndex]);
  if (categoryIndex == 0) {
    y_start = start[1];
  }
  var end = api.coord([api.value(2), categoryIndex]);
  var width = end[0] - start[0];
  var x = start[0];
  var y;
  if (api.value(5) == "gene_feature") {
    y = y_start - height / 2 - api.value(3);
    points = [
      [x, y],
      [x + width - width / 100, y],
      [x + width, y - height / 2],
      [x + width - width / 100, y - height],
      [x, y - height],
    ];
    shape = echarts.graphic.clipPointsByRect(points, {
      x: params.coordSys.x,
      y: params.coordSys.y,
      width: params.coordSys.width,
      height: params.coordSys.height,
    });
    style = api.style();
    textContent = {
      type: "text",
      invisible: false,
      style: {
        text: gene_feature[categoryIndex].name,
        fill: gene_feature[categoryIndex].itemStyle.color,
        fontStyle: "normal",
        fontSize: 10,
        fontWeight: "bolder",
      },
    };
    textConfig = {
      position: "top",
      distance: 15,
      rotation: 0.7,
    };
  } else if (api.value(5) == "amplicon_feature") {
    y = y_start + height / 2 - api.value(3);
    points = [
      [x, y],
      [x + width, y],
      [x + width, y - height],
      [x, y - height],
    ];
    shape = echarts.graphic.clipPointsByRect(points, {
      x: params.coordSys.x,
      y: params.coordSys.y,
      width: params.coordSys.width,
      height: params.coordSys.height,
    });
    style = api.style();
    textContent = {};
    textConfig = {};
  }
  return {
    type: "polygon",
    shape: {
      points: shape,
    },
    style: style,
    textContent: textContent,
    textConfig: textConfig,
  };
}

function updateGeneFeatures(params, api) {
  var categoryIndex = api.value(0);
  var points;
  var shape;
  var style;
  var textConfig;
  var textContent;
  var height = gene_features_properties["rec_items_height"];
  var start = api.coord([api.value(1), categoryIndex]);
  if (categoryIndex == 0) {
    y_start = start[1];
  }
  var end = api.coord([api.value(2), categoryIndex]);
  var width = end[0] - start[0];
  var x = start[0];
  var y;
  if (api.value(5) == "gene_feature") {
    y = y_start - height / 2 - api.value(3);
    points = [
      [x, y],
      [x + width - width / 100, y],
      [x + width, y - height / 2],
      [x + width - width / 100, y - height],
      [x, y - height],
    ];
    shape = echarts.graphic.clipPointsByRect(points, {
      x: params.coordSys.x,
      y: params.coordSys.y,
      width: params.coordSys.width,
      height: params.coordSys.height,
    });
    style = api.style();
    textContent = {
      type: "text",
      invisible: true,
    };
    textConfig = {};
  } else if (api.value(5) == "amplicon_feature") {
    y = y_start + height / 2 - api.value(3);
    points = [
      [x, y],
      [x + width, y],
      [x + width, y - height],
      [x, y - height],
    ];
    shape = echarts.graphic.clipPointsByRect(points, {
      x: params.coordSys.x,
      y: params.coordSys.y,
      width: params.coordSys.width,
      height: params.coordSys.height,
    });
    style = api.style();
    textContent = {};
    textConfig = {};
  }
  return {
    type: "polygon",
    shape: {
      points: shape,
    },
    style: style,
    textContent: textContent,
    textConfig: textConfig,
  };
}

function getGeneFeatureSeries(index) {
  var feature_series = [];
  feature_series.push({
    type: "custom",
    xAxisIndex: index,
    yAxisIndex: index,
    renderItem: renderGeneFeatures,
    data: gene_feature,
    tooltip: {
      trigger: "item",
      enterable: true,
      appendToBody: true,
      renderMode: "html",
      borderRadius: 6,
      borderWidth: 2,
      showContent: "true",
      textStyle: {
        fontSize: 15,
        fontWeight: "bolder",
      },
      formatter: function (params) {
        var output = "";
        output +=
          params.name +
          "<br/>" +
          "Start pos: " +
          params.value[1].toLocaleString() +
          "<br/>" +
          "End pos: " +
          params.value[2].toLocaleString() +
          "<br/>" +
          "Length: " +
          (params.value[2] - params.value[1] + 1).toLocaleString();
        return output;
      },
    },
  });
  return feature_series;
}
///////////////////// End of Gene Feature/////////////////////

function getXAxes(samples, ref_len) {
  var axes = [];
  for (var [i, sample] of samples.entries()) {
    axes.push({
      type: "value",
      gridIndex: i,
      min: 1,
      max: ref_len,
      axisLabel: {
        interval: "auto",
      },
    });
  }
  axes.push({
    type: "value",
    gridIndex: samples.length,
    min: 1,
    max: ref_len,
    axisLabel: {
      interval: "auto",
    },
  });
  return axes;
}

function getYAxes(samples, scaletype, ymax) {
  var axes = [];
  for (var [i, sample] of samples.entries()) {
    axes.push({
      type: scaletype,
      gridIndex: i,
      name: sample,
      nameTextStyle: {
        fontStyle: "normal",
        fontWeight: "bolder",
      },
      nameLocation: "end",
      min: 1,
      max: ymax,
      minorSplitLine: {
        show: true,
      },
    });
  }
  axes.push({
    max: gene_features_properties["max_grid_height"],
    gridIndex: samples.length,
    show: false,
  });
  return axes;
}

function getDatasets(depths, depths_pool1, depths_pool2, positions) {
  var datasets = [];
  for (var [i, depthArray] of depths.entries()) {
    datasets.push({
      dimensions: [
        { name: "depth", type: "float" },
        { name: "depth_pool1", type: "float" },
        { name: "depth_pool2", type: "float" },
        { name: "position", type: "int" },
      ],
      source: {
        position: positions,
        depth: depthArray,
        depth_pool1: depths_pool1[i],
        depth_pool2: depths_pool2[i],
      },
    });
  }
  return datasets;
}

function getDepthSeries(samples) {
  var series = [];
  for (var [i, sample] of samples.entries()) {
    series.push(
      {
        type: "line",
        xAxisIndex: i,
        yAxisIndex: i,
        areaStyle: {
          color: "#666",
        },
        encode: {
          x: "position",
          y: "depth",
        },
        symbol: "none",
        datasetIndex: i,
        lineStyle: {
          color: "#666",
          opacity: 0,
        },
        large: true,
      },
      {
        type: "bar",
        name: "Pool1",
        xAxisIndex: i,
        yAxisIndex: i,
        itemStyle: {
          color: "violet",
        },
        encode: {
          x: "position",
          y: "depth_pool1",
        },
        barWidth: "100%",
        emphasis: {
          focus: "self",
        },
        blur: {},
        datasetIndex: i,
        large: true,
      },
      {
        type: "bar",
        name: "Pool2",
        xAxisIndex: i,
        yAxisIndex: i,
        itemStyle: {
          color: "skyblue",
        },
        encode: {
          x: "position",
          y: "depth_pool2",
        },
        barWidth: "100%",
        emphasis: {
          focus: "self",
        },
        blur: {},
        datasetIndex: i,
        large: true,
      }
    );
  }
  return series;
}

function getVariantsSeries(variants, depths) {
  var series = [];
  for (var [i, varMap] of variants.entries()) {
    (function (i, varMap) {
      series.push({
        type: "bar",
        xAxisIndex: i,
        yAxisIndex: i,
        data: Object.keys(varMap).map((x) => [parseInt(x), depths[i][x]]),
        barWidth: 2,
        itemStyle: {
          color: function (params) {
            var pos = params.data[0];
            var nt = variants[i][pos];
            if (ntColor.hasOwnProperty(nt[0][0])) {
              return ntColor[nt[0][0]];
            } //else {
            //console.log(variants.length, pos, i)
            //}
            return "#333";
          },
        },
      });
    })(i, varMap);
  }
  return series;
}

function getGrids(samples) {
  var n = samples.length + 1;
  var last_height;
  var last_top;
  var grids = Object.keys(samples).map(function (sample) {
    last_height = (1 / n) * 100 - 6;
    if (n == 2) {
      // Only 1 sample (1 sample + gene feature plot)
      last_height = 70;
      return {
        show: true,
        height: "70%", // plot display in nearly full scale
      };
    }
    return {
      show: true,
      height: (1 / n) * 100 - 6 + "%",
    };
  });
  grids.forEach(function (grid, idx) {
    //var padTop = idx === 1 ? 5 : 3
    var padTop = 4;
    last_top = (idx / n) * 100 + padTop;
    grid.top = (idx / n) * 100 + padTop + "%";
    grid.left = "5%";
    grid.right = "5%";
  });
  grids.push({
    show: true,
    height: gene_features_properties["grid_height"],
    top: last_height + last_top + 3 + "%",
    left: "5%",
    right: "5%",
  });
  return grids;
}

function meanCoverage(depths, start, end, gridIndex) {
  var total = 0;
  if (start < end) {
    for (var i = start - 1; i <= end - 1; i++) {
      if (depths[gridIndex][i] === 1e-6) {
        continue;
      }
      total += depths[gridIndex][i];
    }
    return total / (end - start + 1);
  } else if (start == end) {
    if (depths[gridIndex][start] === 1e-6) {
      return 0;
    } else {
      return depths[gridIndex][start];
    }
  }
}

function genomeCoverage(depths, start, end, gridIndex, low) {
  var total = 0;
  for (var i = start - 1; i <= end - 1; i++) {
    if (depths[gridIndex][i] >= low) {
      total += 1;
    }
  }
  return (total / (end - start + 1)) * 100;
}

function median(arr) {
  arr.sort(function (a, b) {
    return a - b;
  });
  var half = Math.floor(arr.length / 2);
  if (arr.length % 2) {
    if (arr[half] === 1e-6) {
      return 0;
    } else {
      return arr[half];
    }
  } else if (arr[half - 1] != 1e-6 && arr[half] != 1e-6)
    return (arr[half - 1] + arr[half]) / 2.0;
  else {
    if (arr[half - 1] === 1e-6) {
      return arr[half] / 2.0;
    } else if (arr[half] === 1e-6) {
      return arr[half - 1] / 2.0;
    }
  }
}

function medianCoverage(depths, start, end, gridIndex) {
  var sub_array = depths[gridIndex].slice(start - 1, end);
  return median(sub_array);
}

function getTooltips(samples, depths, variants) {
  return [
    {
      trigger: "axis",
      enterable: true,
      appendToBody: true,
      renderMode: "html",
      showContent: true,
      formatter: function (params) {
        var output = "";
        var param = params[0];
        var i = param.axisIndex;
        var position = param.data[3];
        if (i < samples.length) {
          var sample = samples[i];
          var depth = param.data[0];
          var amplicon_depth_pool1 = param.data[1];
          var amplicon_depth_pool2 = param.data[2];
          var start_pos, end_pos;
          start_pos = Math.floor(chart.getOption().dataZoom[0].startValue);
          end_pos = Math.floor(chart.getOption().dataZoom[0].endValue);
          var mean_cov;
          mean_cov = meanCoverage(depths, start_pos, end_pos, i).toFixed(2);
          var median_cov;
          median_cov = medianCoverage(depths, start_pos, end_pos, i).toFixed(2);
          var genome_cov;
          genome_cov = genomeCoverage(
            depths,
            start_pos,
            end_pos,
            i,
            10
          ).toFixed(2);
          output +=
            "<b>" +
            sample +
            "</b><br/>" +
            "Start pos: " +
            start_pos.toLocaleString() +
            "<br/>" +
            "End pos: " +
            end_pos.toLocaleString() +
            "<br/>";
          output += "Genome Mean Coverage: " + mean_cov + "X" + "<br/>";
          output += "Genome Median Coverage: " + median_cov + "X" + "<br/>";
          output += "Genome Coverage ( >= 10x): " + genome_cov + "%" + "<br/>";
          output +=
            "Genome Position: " +
            position.toLocaleString() +
            "<br/> Genome Coverage Depth: " +
            depth.toLocaleString() +
            "<br/>";
          output +=
            "Pool1 Coverage Depth: " +
            amplicon_depth_pool1.toLocaleString() +
            "<br/>";
          output +=
            "Pool2 Coverage Depth: " +
            amplicon_depth_pool2.toLocaleString() +
            "<br/>";
          if (variants[i].hasOwnProperty(position)) {
            output +=
              "Ref: " +
              window.ref_seq.substring(
                position - 1,
                position - 1 + variants[i][position].length
              ) +
              "<br/>";
            output = output + "Variant: " + variants[i][position];
          } else {
            output += "Ref: " + window.ref_seq[position - 1] + "<br/>";
          }
          return output;
        }
        return output;
      },
    },
  ];
}

function getOption() {
  var samples = [];
  var depths = [];
  var depths_pool1 = [];
  var depths_pool2 = [];
  var variants = [];

  for (const [key, entries] of Object.entries(window.samples)) {
    if (key < default_num_chart) {
      samples.push(entries);
      depths.push(window.depths[entries][0]);
      depths_pool1.push(window.depths[entries][1]);
      depths_pool2.push(window.depths[entries][2]);
      variants.push(window.variants[entries]);
    }
  }
  grid_length = samples.length;
  var options = {
    title: {},
    dataset: getDatasets(depths, depths_pool1, depths_pool2, positions),
    xAxis: getXAxes(samples, ref_len),
    yAxis: getYAxes(samples, "log", 100000),
    legend: {
      data: ["Pool1", "Pool2"],
    },
    series: [
      ...getDepthSeries(samples),
      ...getVariantsSeries(variants, depths),
      ...getGeneFeatureSeries(grid_length),
    ],
    tooltip: getTooltips(samples, depths, variants),
    toolbox: {
      show: "true",
      feature: {
        restore: {},
        saveAsImage: {
          name: "Coverage_Plot",
        },
      },
    },
    dataZoom: [
      {
        type: "inside",
        filterMode: "none",
        xAxisIndex: [...Array(grid_length + 1).keys()],
      },
      {
        show: true,
        filterMode: "none",
        xAxisIndex: [...Array(grid_length + 1).keys()],
        type: "slider",
      },
    ],
    grid: getGrids(samples),
  };
  selectDefaultSamples(samples);
  return options;
}

function updateOption(samples) {
  var depths = [];
  var depths_pool1 = [];
  var depths_pool2 = [];
  var variants = [];

  for (const selected_samples of samples) {
    depths.push(window.depths[selected_samples][0]);
    depths_pool1.push(window.depths[selected_samples][1]);
    depths_pool2.push(window.depths[selected_samples][2]);
    variants.push(window.variants[selected_samples]);
  }
  grid_length = samples.length;
  var options = {
    title: {},
    dataset: getDatasets(depths, depths_pool1, depths_pool2, positions),
    xAxis: getXAxes(samples, ref_len),
    yAxis: getYAxes(samples, "log", 100000),
    legend: {
      data: ["Pool1", "Pool2"],
    },
    series: [
      ...getDepthSeries(samples),
      ...getVariantsSeries(variants, depths),
      ...getGeneFeatureSeries(grid_length),
    ],
    tooltip: getTooltips(samples, depths, variants),
    toolbox: {
      show: "true",
      feature: {
        restore: {},
        saveAsImage: {
          name: "Coverage_Plot",
        },
      },
    },
    dataZoom: [
      {
        type: "inside",
        filterMode: "none",
        xAxisIndex: [...Array(grid_length + 1).keys()],
      },
      {
        show: true,
        filterMode: "none",
        xAxisIndex: [...Array(grid_length + 1).keys()],
        type: "slider",
      },
    ],
    grid: getGrids(samples),
  };
  chart.setOption((option = options), (notMerge = true));
  updateChartOptionsMenu();
}

chart.setOption((option = getOption()));
