function updateChartOptionsMenu() {
  var grid_option = chart.getOption().grid;
  var height1 = _.replace(grid_option[0].height, "%", "");
  var height2 = _.replace(grid_option[grid_option.length - 1].height, "%", "");
  var left = _.replace(grid_option[0].left, "%", "");
  var right = _.replace(grid_option[0].right, "%", "");
  var top = _.replace(grid_option[0].top, "%", "");
  // Reset Chart Properties and Gene Features Menu
  document.getElementById("chartheightinput").value = _.toInteger(height1);
  document.getElementById("chartheightoutput").value = _.toInteger(height1) + "%";
  document.getElementById("chartleftinput").value = _.toInteger(left);
  document.getElementById("chartleftoutput").value = _.toInteger(left) + "%";
  document.getElementById("chartrightinput").value = _.toInteger(right);
  document.getElementById("chartrightoutput").value = _.toInteger(right) + "%";
  document.getElementById("charttopinput").value = _.toInteger(top);
  document.getElementById("charttopoutput").value = _.toInteger(top) + "%";
  document.getElementById("genefeatureheightinput").value = _.toInteger(height2);;
  document.getElementById("genefeatureheightoutput").value = _.toInteger(height2); + "%";
  // Set Axis to Log scale
  document.getElementById("scale").value = "log"
  document.getElementById("ymax").value = 100000
  // Reset Data Zoom
  resetDataZoom()
}

function chartDispatchAction() {
  chart.on("click", function (params) {
    document.getElementById("start_pos").value = params.value[1];
    document.getElementById("end_pos").value = params.value[2];
    chart.dispatchAction({
      type: "dataZoom",
      startValue: params.value[1],
      endValue: params.value[2],
    });
  });

  chart.on("dblclick", function (params) {
    document.getElementById("start_pos").value = 1;
    document.getElementById("end_pos").value = ref_len;
    chart.dispatchAction({
      type: "dataZoom",
      start: 0,
      end: 100,
    });
  });
}