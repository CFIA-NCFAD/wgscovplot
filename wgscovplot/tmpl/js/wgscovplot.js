/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are reserved
 * Users's settings are respected
 * @param {Array<string>} samples - An array of samples name
 */
function updateCoverageChartOption(samples) {
    var depths = [];
    var variants = [];
    var scaleType;
    var max;
    var leftMargin;
    var rightMargin;
    samples.forEach(sample => {
        depths.push(window.depths[sample]);
        variants.push(window.variants[sample]);
    })
    /* Need to handle the case when user removes all samples on the chart and add new samples to be displayed
    Variables lenCheck1, lenCheck2 and closures lastYAxisScale, lastYAxisMax are used to keep track the scale settings of a last sample on the chart before it was removed
    lenCheck1, lenCheck1 depends on whether the last subplot is used for plotting gene/amplicon features or not
     */
    var chartOption = chart.getOption();
    var lenCheck1 = (amplicon === 'True' || geneFeature === 'True') ? 2 : 1;
    var lenCheck2 = (amplicon === 'True' || geneFeature === 'True') ? 1 : 0;
    if (chartOption.yAxis.length === lenCheck1){
        lastYAxisScale = chartOption.yAxis[0].type;
        lastYAxisMax = chartOption.yAxis[0].max;
    }
    if (chartOption.yAxis.length === lenCheck2){
        scaleType = lastYAxisScale
        max = lastYAxisMax
    }else{
        // Get Current scale type and yAxis max to set it back
        scaleType = chartOption.yAxis[0].type;
        max = chartOption.yAxis[0].max;
    }
    // Left-right margin
    if (chartOption.grid.length === 0){ // There is no subplot and gene/amplicon feature
        leftMargin = 8
        rightMargin = 8
        document.getElementById("chart-left-input").value = parseInt(leftMargin);
        document.getElementById("chart-left-output").value = parseInt(leftMargin) + "%";
        document.getElementById("chart-right-input").value = parseInt(rightMargin);
        document.getElementById("chart-right-output").value = parseInt(rightMargin) + "%";
    }
    else{
        leftMargin = parseInt(chartOption.grid[0].left.replace("%", ""));
        rightMargin = parseInt(chartOption.grid[0].right.replace("%", ""));
    }
    // Get Current Data Zoom
    var zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
    var zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
    // The current chart is not disposed so notMerge must be set true
    chart.setOption(option = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, positions,
        maxDepth, samples, depths, variants, geneFeature, amplicon), notMerge = true);
    setScale(scaleType, max);
    setDataZoom(zoomStart, zoomEnd);
    updateChartLeftMargin(leftMargin);
    updateChartRightMargin(rightMargin);
    updateControlMenu();
}

/**
 * When the chart is initialized, the first 3 samples are plotted
 * @param {Array<string>} samples - An array of samples name
 */
function setDefaultSamples(samples) {
    // Set default samples display
    let $selectedsamples = $("#selectedsamples");
    $selectedsamples.select2();
    $selectedsamples.val(samples);
    $selectedsamples.trigger('change')
}

/**
 * Get the list of current samples on the control menu
 * @returns {Array<string>>} An array of samples name
 */
function getCurrentSamples(chartOption) {
    var samples = [];
    // If the chart is not initalizaed yet, get 3 first samples from window.samples
    if (chartOption === undefined || chartOption === null){
        samples = window.samples.slice(0, 3)
    } else{
        var selectData = $("#selectedsamples").select2("data");
        for (var [key, entries] of selectData.entries()) {
            samples.push(selectData[key].text);
        }
    }
    return samples;
}

/**
 * Initialize the events handler for chart
 */
function initWgscovplotEvent(){
    $(document).ready(function () {
        /**
         * Toggle dark mode, the entire chart will be disposed and re-initialized
         * However, the old settings of charts are set back (users's settings are respected)
         */
        $("#toggle-darkmode").change(function () {
            initWgscovplotRenderEnv();
        });

        /**
         * Toggle to Amplicon Depth Label
         */
        $("#toggle-amplicon-depthlabel").change(function () {
            var seriesOption = chart.getOption().series;
            var isChecked = $(this).prop("checked");
            seriesOption.forEach(element => {
                if (element.type === 'custom') {
                    element.label.show = isChecked
                }
            });
            chart.setOption({series: [...seriesOption]});
        });

        /**
         * Jquery actions to make the list of samples is not forced in alphabetical order
         */
        $("#selectedsamples").select2({
            tags: true,
        });

        $("#selectedsamples").on("select2:select", function (evt) {
            var element = evt.params.data.element;
            var $element = $(element);
            $element.detach();
            $(this).append($element);
            $(this).trigger("change");
        });

        /**
         * Update chart options when the number of selected samples changed
         * The chart options such as y Axis scale, yMax, DataZoom are reserved (users's settings are respected)
         */
        $("#selectedsamples").on("change", function () {
            updateCoverageChartOption(getCurrentSamples(chart.getOption()));
        });

        /**
         * Toggle to show gene label or not
         */
        $("#toggle-genelabel").change(function () {
            var seriesOption = chart.getOption().series;
            showGeneLabel = $(this).prop("checked");
            seriesOption[seriesOption.length - 1]["renderItem"] = wgscovplot.getGeneFeatureRenderer(showGeneLabel, geneFeatureAmpliconData); // Re-update Gene Feature Chart Only
            chart.setOption({series: [...seriesOption]});
        });

        /**
         * Toggle tooltip for coverage chart, this action does not affect gene/amplicon feature chart
         */
        $("#toggle-tooltip").change(function () {
            var isChecked = $(this).prop("checked");
            chart.setOption({
                tooltip: {showContent: isChecked},
            });
        });

        /**
         * Toggle slider zoom
         */
        $("#toggle-slider").change(function () {
            var isChecked = $(this).prop("checked");
            var numChart = chart.getOption().grid.length;
            chart.setOption({
                dataZoom: [
                    {
                        type: "inside",
                        filterMode: "none",
                        xAxisIndex: [...Array(numChart).keys()],
                    },
                    {
                        type: "slider",
                        show: isChecked,
                        filterMode: "none",
                        xAxisIndex: isChecked ? [...Array(numChart).keys()] : null,
                    },
                ],
            })
        });
    });
}

/**
 * Initialize the enviroment for the chart (sgv/canvas or dark/white mode)
 * The entire chart will be disposed and re-initialized
 * However, the old settings of charts are reserved (users's settings are respected)
 */
function initWgscovplotRenderEnv() {
    var chartOption = chart.getOption();
    var plotSamples = getCurrentSamples(chartOption);
    var plotDepths = [];
    var plotVariants = [];
    plotSamples.forEach(sample => {
        plotDepths.push(window.depths[sample]);
        plotVariants.push(window.variants[sample]);
    });
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamples(plotSamples);
        chart.setOption(option = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, positions,
            maxDepth, plotSamples, plotDepths, plotVariants, geneFeature, amplicon));
    } else {
        var renderEnv = document.getElementById("render-env").value;
        var isChecked = document.getElementById("toggle-darkmode").checked;
        var mode = isChecked ? "dark" : "white";
        var gridOption = chart.getOption().grid;
        var scaleType = chart.getOption().yAxis[0].type;
        var max = chart.getOption().yAxis[0].max;
        var dataZoomOption = chart.getOption().dataZoom;
        wgscovplot.echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        chart = wgscovplot.echarts.init($chart, mode, {renderer: renderEnv});
        chart.setOption(option = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, positions,
            maxDepth, plotSamples, plotDepths, plotVariants, geneFeature, amplicon));
        chart.setOption({grid: gridOption, dataZoom: dataZoomOption});
        setScale(scaleType, max);
    }
    onChartDataZoomActions();
}

/**
 * Adjust subplot height
 * @param {number} val - Subplots height percent value
 */
function updateSubPlotHeight(val) {
    document.getElementById("chart-height-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    var len = (amplicon === 'True' || geneFeature === 'True') ? gridOption.length - 1 : gridOption.length
    for (var i = 0; i < len; i++) {
        gridOption[i]["height"] = val + "%";
        if (i > 0) {
            gridOption[i]["top"] =
                parseInt(gridOption[i - 1]["top"].replace("%", "")) +
                parseInt(gridOption[i - 1]["height"].replace("%", "")) +
                4 +
                "%";
        }
    };
    chart.setOption({grid: gridOption});
}

/**
 * Adjust left margin of chart
 * @param {number} val - Left margin percent value
 */
function updateChartLeftMargin(val) {
    document.getElementById("chart-left-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    gridOption.forEach(element => {
        element.left = val + "%";
    });
    chart.setOption({grid: gridOption});
}

/**
 * Adjust right margin of chart
 * @param {number} val - Right margin percent value
 */
function updateChartRightMargin(val) {
    document.getElementById("chart-right-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    gridOption.forEach(element => {
        element.right = val + "%";
    });
    chart.setOption({grid: gridOption});
}

/**
 * Adjust top margin of chart
 * @param {number} val - Subplots top margin percent value
 */
function updateSubPlotTopMargin(val) {
    document.getElementById("chart-top-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    for (var i = 0; i < gridOption.length; i++) {
        if (i === 0) {
            gridOption[i]["top"] = val + "%";
        } else {
            gridOption[i]["top"] =
                parseInt(gridOption[i - 1]["top"].replace("%", "")) +
                parseInt(gridOption[i - 1]["height"].replace("%", "")) +
                parseInt(val) +
                "%";
        }
    }
    chart.setOption({grid: gridOption});
}

/**
 * Adjust the height of gene/amplicon feature charts
 * @param {number} val - Subplots gene feature height percent value
 */
function updateGeneFeatureHeight(val) {
    document.getElementById("genefeature-height-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    gridOption[gridOption.length - 1]["height"] = val + "%";
    chart.setOption({grid: gridOption});
}

/**
 * Set scale for y Axis
 * @param {string} scaleType - Either value or log
 * @param {number} max - Max value is set for y Axis
 */
function setScale(scaleType, max) {
    var scale;
    var yMax;
    if (scaleType === null && max === null){
        scale = document.getElementById("scale").value;
        yMax = document.getElementById("ymax").value;
    }
    else{
        scale = scaleType;
        yMax = max;
    }
    var yAxisOption = chart.getOption().yAxis;
    var len = (amplicon === 'True' || geneFeature === 'True') ? yAxisOption.length - 1 : yAxisOption.length
    if (scale === "value") {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scale;
                element.min = 0;
                element.max = yMax;
            }
        });
    } else {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scale;
                element.min = 1;
                element.max = yMax;
            }
        });
    }
    chart.setOption({yAxis: yAxisOption});
    document.getElementById("ymax").value = yMax; // update control menu
}

/**
 * Set yMax for Y Axis
 */
function setYMax() {
    var yMax = document.getElementById("ymax").value;
    var yAxisOption = chart.getOption().yAxis;
    var len = (amplicon === 'True' || geneFeature === 'True') ? yAxisOption.length - 1 : yAxisOption.length
    yAxisOption.forEach(element => {
        if (element.gridIndex < len) {
            element.max = yMax;
        }
    });
    chart.setOption({yAxis: yAxisOption});
}

/**
 * Set zoom view for the chart
 * @param {number} zoomStart - Start view point
 * @param {number} zoomEnd - End view point
 */
function setDataZoom(zoomStart, zoomEnd){
    var start;
    var end;
    if (zoomStart === null && zoomEnd === null){
        start = document.getElementById("start-pos").value;
        end = document.getElementById("end-pos").value;
    }
    else{
        start = zoomStart;
        end = zoomEnd;
    }
    // Update Control menu
    document.getElementById("start-pos").value = start;
    document.getElementById("end-pos").value = end;
    // Dispatch Action to chart
    chart.dispatchAction({
        type: "dataZoom",
        startValue: start,
        endValue: end,
    });
}

/**
 * Dispatch click/dbclick actions for the whole chart
 */
function onChartDataZoomActions(){
    chart.on("click", function(params){
        document.getElementById("start-pos").value = params.value.start;
        document.getElementById("end-pos").value = params.value.end;
        setDataZoom(params.value.start, params.value.end);
    });
    chart.on("dblclick", function(params){
        document.getElementById("start-pos").value = 1;
        document.getElementById("end-pos").value = refSeqLength;
        setDataZoom(1, refSeqLength);
    });
}

/**
 * The Control Menu is updated when the number of selected samples changes
 * Menu is updated to reflect chart properties such as subplot height/top margin
 */
function updateControlMenu() {
    var gridOption = chart.getOption().grid;
    if (gridOption.length > 0) {
        var height1 = gridOption[0].height.replace("%", "");
        var top = gridOption[0].top.replace("%", "");
        document.getElementById("chart-height-input").value = parseInt(height1);
        document.getElementById("chart-height-output").value = parseInt(height1) + "%";
        document.getElementById("chart-top-input").value = parseInt(top);
        document.getElementById("chart-top-output").value = parseInt(top) + "%";
        if (amplicon === 'True' || geneFeature === 'True') {
            var height2 = gridOption[gridOption.length - 1].height.replace("%", "");
            document.getElementById("genefeature-height-input").value = parseInt(height2);
            document.getElementById("genefeature-height-output").value = parseInt(height2) + "%";
        }
    }
}

