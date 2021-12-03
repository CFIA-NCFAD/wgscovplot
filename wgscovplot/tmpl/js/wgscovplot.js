/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are reserved
 * Users's settings are respected by keeping old settings and set it back.
 * @param {Array<string>} samples - An array of samples name
 */
function updateCoverageChartOption(samples) {
    var depths = [];
    var variants = [];
    var scaleType;
    var yAxisMax;
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
        yAxisMax = lastYAxisMax
    }else{
        // Get Current scale type and yAxis max to set it back
        scaleType = chartOption.yAxis[0].type;
        yAxisMax = chartOption.yAxis[0].max;
    }
    // Left-right margin
    if (chartOption.grid.length === 0){ // There is no subplot and gene/amplicon feature
        leftMargin = "8%";
        rightMargin = "8%";
        document.getElementById("chart-left-input").value = parseInt(leftMargin);
        document.getElementById("chart-left-output").value = parseInt(leftMargin) + "%";
        document.getElementById("chart-right-input").value = parseInt(rightMargin);
        document.getElementById("chart-right-output").value = parseInt(rightMargin) + "%";
    }
    else{
        leftMargin = chartOption.grid[0].left;
        rightMargin = chartOption.grid[0].right;
    }

    // Get Current Data Zoom
    var oldDataZoom = chartOption.dataZoom;

    // get Coverage Chart Option with new data
    var updateOption = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, window.refSeq,
        yAxisMax, samples, depths, variants, geneFeature, amplicon);

    // Update tooltip
    var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
    var triggerOnType;
    if (isTooltipEnable){
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    }else{
        triggerOnType ="none";
    }
    var isVariantComparison = document.getElementById("toggle-variant-comparison").checked;
    var isVariantSitesOnly = document.getElementById("toggle-tooltip-variant-sites").checked;
    updateOption.tooltip = wgscovplot.getTooltips(samples, depths, variants, window.refSeq,
                triggerOnType=triggerOnType, isVariantSitesOnly=isVariantSitesOnly,
                isVariantComparison=isVariantComparison)
    var isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
    updateOption.tooltip[0]["position"] = tooltipPosition(isFixTooltipPostion);

    // Update grid
    updateOption.grid.forEach(element => {
        element.left = leftMargin;
        element.right = rightMargin;
    })

    // Update datzoom
    updateOption.dataZoom = oldDataZoom;
    updateOption.dataZoom.forEach(element => {
        element.xAxisIndex = [...Array(updateOption.grid.length).keys()];
    })

    //update scale and yAxis max
    updateOption.yAxis = updateYAxisOption(updateOption.yAxis, scaleType, yAxisMax);

    //set chart option
    chart.setOption(option = updateOption, notMerge = true);

    //update control menu
    updateControlMenu();
}

/**
 * Update scale type and max for Y Axis
 * @param {Dict[]} yAxisOption - Options of Yaxis need to be updated
 * @param {string} scaleType - Either log or value
 * @param {number} yAxisMax - Max value is set for Y Axis
 * @returns {Dict[]} Returns the updated options (Scale type or ymax) for yAxis
 */
function updateYAxisOption(yAxisOption, scaleType, yAxisMax){
    var len = (amplicon === 'True' || geneFeature === 'True') ? yAxisOption.length - 1 : yAxisOption.length
    if (scaleType === "value") {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 0;
                element.max = yAxisMax;
            }
        });
    } else {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 1;
                element.max = yAxisMax;
            }
        });
    }
    return yAxisOption;
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
         * Update chart options when adding/romoving samples
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
         * Toggle tooltip for coverage chart
         */
        $("#toggle-tooltip").change(function () {
            var isChecked = $(this).prop("checked");
            var triggerType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
            chart.setOption({
                tooltip: {triggerOn: isChecked ? triggerType: "none"},
            });
        });

        /**
         * Toggle turn tooltip display for variant sites only
         */
        $("#toggle-tooltip-variant-sites").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var samples = getCurrentSamples(chartOption);
            var depths = [];
            var variants = [];
            samples.forEach(sample => {
                depths.push(window.depths[sample]);
                variants.push(window.variants[sample]);
            });
            var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
            var triggerOnType;
            if (isTooltipEnable){
                triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
            }else{
                triggerOnType ="none";
            }
            var isVariantComparison = document.getElementById("toggle-variant-comparison").checked;
            var tooltipOption = wgscovplot.getTooltips(samples, depths, variants, window.refSeq,
                triggerOnType=triggerOnType, isVariantSitesOnly=isChecked,
                isVariantComparison=isVariantComparison);
            var isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
            tooltipOption[0]["position"] = tooltipPosition(isFixTooltipPostion);
            chart.setOption({tooltip: tooltipOption});
        });


        $("#toggle-variant-comparison").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var samples = getCurrentSamples(chartOption);
            var depths = [];
            var variants = [];
            samples.forEach(sample => {
                depths.push(window.depths[sample]);
                variants.push(window.variants[sample]);
            });
            var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
            var triggerOnType;
            if (isTooltipEnable){
                triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
            }else{
                triggerOnType ="none";
            }
            var isVariantSitesOnly = document.getElementById("toggle-tooltip-variant-sites").checked;
            var tooltipOption = wgscovplot.getTooltips(samples, depths, variants, window.refSeq,
                triggerOnType=triggerOnType, isVariantSitesOnly=isVariantSitesOnly,
                isVariantComparison=isChecked);
            var isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
            tooltipOption[0]["position"] = tooltipPosition(isFixTooltipPostion);
            chart.setOption({tooltip: tooltipOption});
        });

        /**
         * Fix tooltip position on the charts
         * If it is fixed: tooltip will be on the right if mouse hovering on the left and vice versa
         * Default tooltip follows cursor
         */
        $("#fix-tooltip-postion").change(function () {
            var isChecked = $(this).prop("checked");
            var tooltipOption  = chart.getOption().tooltip;
            tooltipOption[0]["position"] = tooltipPosition(isChecked)
            chart.setOption({tooltip: tooltipOption});
        });

        /**
         * Toogle tooltip trigger on click
         */
        $("#toggle-tooltip-trigger-click").change(function () {
            var isChecked = $(this).prop("checked");
            var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
            if (isTooltipEnable){
                chart.setOption({
                    tooltip: {triggerOn: isChecked ? "click" : "mousemove"}
                });
            }
            else{
                chart.setOption({
                    tooltip: {triggerOn: "none"}
                });
            }
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
        chart.setOption(option = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, window.refSeq,
            yAxisMax, plotSamples, plotDepths, plotVariants, geneFeature, amplicon,
            triggerOnType= "mousemove", isVariantSitesOnly = false, isVariantComparison=true));
    } else {
        var renderEnv = document.getElementById("render-env").value;
        var isChecked = document.getElementById("toggle-darkmode").checked;
        var mode = isChecked ? "dark" : "white";
        var gridOption = chart.getOption().grid;
        var scaleType = chart.getOption().yAxis[0].type;
        var yAxisMax = chart.getOption().yAxis[0].max;
        var dataZoomOption = chart.getOption().dataZoom;
        var tooltipOption = chart.getOption().tooltip;
        wgscovplot.echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        chart = wgscovplot.echarts.init($chart, mode, {renderer: renderEnv});
        var chartOption = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, window.refSeq,
            yAxisMax, plotSamples, plotDepths, plotVariants, geneFeature, amplicon)
        // Keep grid option
        chartOption.grid = gridOption;
        // Keep data zoom option
        chartOption.dataZoom = dataZoomOption;
        // Keep yAxis option
        chartOption.yAxis = updateYAxisOption(chartOption.yAxis, scaleType, yAxisMax);
        // Keep tooltip
        chartOption.tooltip = tooltipOption;
        //set chart option
        chart.setOption(option = chartOption);
    }
    onChartDataZoomActions();
}

/**
 * Functions to return tooltip position
 * @param {boolean} isChecked
 * @returns Returns the tooltip position
 */
function tooltipPosition(isChecked){
    if (isChecked){
        return function (pos, params, dom, rect, size) {
            var obj = {top: 5};
            obj[['left', 'right'][+(pos[0] < size.viewSize[0] / 2)]] = 5;
            return obj;
        }
    }
    else{
        return ""; // follow cursor
    }
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
 */
function setScale() {
    var scaleType = document.getElementById("scale").value;
    var yAxisMax = document.getElementById("ymax").value;
    var yAxisOption =  updateYAxisOption(chart.getOption().yAxis, scaleType, yAxisMax)
    chart.setOption({yAxis: yAxisOption});
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

