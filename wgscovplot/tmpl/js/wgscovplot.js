/**
 * Get depths and variants for samples
 * @param {Array<string>} samples - An array of samples name
 * @param {Object} depthsData - Depths data
 * @param {Object} variantsData - Variants data
 * @returns [Array<Array<number>>, Array<Array<Object>>]
 */
function getDepthsVariants(samples, depthsData, variantsData){
    var depths = [];
    var variants = [];
    samples.forEach(sample => {
        if (depthsData[sample] !== undefined && depthsData[sample] !== null)
            depths.push(depthsData[sample]);
        else
            depths.push([]);
        if (variantsData[sample] !== undefined && variantsData[sample] !== null)
            variants.push(variantsData[sample]);
        else
            variants.push({});
    })
    return [depths, variants]
}
/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are reserved
 * Users's settings are respected by keeping old settings and set it back.
 * @param {Array<string>} samples - An array of samples name
 */
function updateCoverageChartOption(samples) {
    var scaleType;
    var yAxisMax;
    const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
    var chartOption = chart.getOption();
    // Reserver Tooltip option
    var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
    var triggerOnType;
    if (isTooltipEnable){
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    }else{
        triggerOnType ="none";
    }
    var isVariantSites = document.getElementById("toggle-tooltip-variant-sites").checked;
    var isNonVarianSites = document.getElementById("toggle-tooltip-non-variant-sites").checked;
    var isVariantComparison = document.getElementById("toggle-variant-comparison").checked;
    var isCoverageStatView = document.getElementById("toggle-coverage-stat").checked;
    var updateOption = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, window.refSeq,
        yAxisMax, samples, depths, variants, geneFeature, amplicon,
        triggerOnType= triggerOnType, isVariantSites=isVariantSites,
        isNonVariantSites=isNonVarianSites, isInfoComparison=isVariantComparison,
        isCovergateStatView=isCoverageStatView);

    // Reserve tooltip in series option
    var seriesOption = updateOption.series;
    seriesOption.forEach(element => {
        if (element.type === 'line'){
            element.tooltip.trigger = isNonVariantSites ? "axis" : "none"
        }
        else if (element.type === 'bar'){
            element.tooltip.trigger = isVariantSites ? "axis" : "none"
        }
    })
    var isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
    updateOption.tooltip[0]["position"] = tooltipPosition(isFixTooltipPostion);

    // Reserve grid option
    updateOption.grid.forEach(element => {
        element.left = $("#chart-left-input").val() + "%";;
        element.right = $("#chart-right-input").val() + "%";;
    })

    // Reserve datzoom
    var oldDataZoom = chartOption.dataZoom;
    updateOption.dataZoom = oldDataZoom;
    updateOption.dataZoom.forEach(element => {
        element.xAxisIndex = [...Array(updateOption.grid.length).keys()];
    })

    //Reserve scale and yAxis max
    scaleType = $("#scale").val();
    yAxisMax = $("#ymax").val();
    updateOption.yAxis = updateYAxisOption(updateOption.yAxis, scaleType, yAxisMax);

    //set chart option
    chart.setOption(option = updateOption, notMerge = true);

    // Update control menu
    updateControlMenu();
}

/**
 * Update scale type and max for Y Axis
 * @param {Object} yAxisOption - Options of Yaxis need to be updated
 * @param {string} scaleType - Either log or value
 * @param {number} yAxisMax - Max value is set for Y Axis
 * @returns {Object} Returns the updated options (Scale type or ymax) for yAxis
 */
function updateYAxisOption(yAxisOption, scaleType, yAxisMax){
    var len = (amplicon || geneFeature) ? yAxisOption.length - 1 : yAxisOption.length
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
 * @returns {Array<string>} An array of samples name
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
         * Toggle turn tooltip display for variant sites only
         */
        $("#toggle-tooltip-variant-sites").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var seriesOption = chartOption.series;
            var samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                isChecked, document.getElementById("toggle-tooltip-non-variant-sites").checked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });

        $("#toggle-tooltip-non-variant-sites").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var seriesOption = chartOption.series;
            var samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked, isChecked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });


        $("#toggle-variant-comparison").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var seriesOption = chartOption.series;
            var samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked,
                document.getElementById("toggle-tooltip-non-variant-sites").checked,
                isChecked, document.getElementById("toggle-coverage-stat").checked);
        });

        $("#toggle-coverage-stat").change(function (){
            var isChecked = $(this).prop("checked");
            var chartOption = chart.getOption();
            var seriesOption = chartOption.series;
            var samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked,
                document.getElementById("toggle-tooltip-non-variant-sites").checked,
                document.getElementById("toggle-variant-comparison").checked, isChecked);
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
    const [plotDepths, plotVariants] = getDepthsVariants(plotSamples, window.depths, window.variants);
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamples(plotSamples);
        chart.setOption(option = wgscovplot.getCoverageChartOption(geneFeatureAmpliconData, ampliconDepthBarData, window.refSeq,
            maxDepth, plotSamples, plotDepths, plotVariants, geneFeature, amplicon));
        variantHeatmap.setOption(option = wgscovplot.getVariantHeatmapOption(window.samples, window.mutations, window.variantMatrix, window.variants));
    } else {
        var renderEnv = document.getElementById("render-env").value;
        var isChecked = document.getElementById("toggle-darkmode").checked;
        var mode = isChecked ? "dark" : "white";
        var gridOption = chart.getOption().grid;
        var seriesOption = chart.getOption().series;
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
        // Keep series
        chartOption.series = seriesOption;
        //set chart option
        chart.setOption(option = chartOption);
    }
    updateControlMenu();
    onChartDataZoomActions();
}

/**
 * Update option for displaying tooltip for variant and non variant sites
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The object of variants data
 * @param {Array<Array<Object>>} seriesOption - Series option of the chart
 * @param {boolean} isVariantSites - Whether to enable tooltip for variant sites or not
 * @param {boolean} isNonVariantSites - Whether to enable tooltip for non variant sites or not
 * @param {boolean} isInfoComparison - Whether to compare variant/ Coverage Stat information across selected samples
 */
function updateTooltipOption(samples, depths, variants, seriesOption,
                              isVariantSites, isNonVariantSites, isInfoComparison, isCovergateStatView){
    var isTooltipEnable = document.getElementById("toggle-tooltip").checked;
    var triggerOnType;
    if (isTooltipEnable){
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    }else{
        triggerOnType ="none";
    }
    seriesOption.forEach(element => {
        if (element.type === 'line'){
            element.tooltip.trigger = isNonVariantSites ? "axis" : "none"
        }
        else if (element.type === 'bar'){
            element.tooltip.trigger = isVariantSites ? "axis" : "none"
        }
    })
    var tooltipOption = wgscovplot.getTooltips(samples, depths, variants, window.refSeq,
    triggerOnType=triggerOnType, isInfoComparison=isInfoComparison, isCovergateStatView=isCovergateStatView);
    var isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
    tooltipOption[0]["position"] = tooltipPosition(isFixTooltipPostion);
    tooltipOption[0].triggerOn = triggerOnType;
    chart.setOption({tooltip: tooltipOption, series: seriesOption});
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
        return "cursor"; // follow cursor
    }
}

/**
 * Adjust Variant Heatmap height
 * @param {number} val - Subplots height percent value
 */
function updateVarMapHeight(val) {
    document.getElementById("varmap-height-output").value = val + "%";
    var gridOption = variantHeatmap.getOption().grid;
    gridOption[0]["height"] = val + "%";
    variantHeatmap.setOption({grid: gridOption});
}

/**
 * Adjust subplot height
 * @param {number} val - Subplots height percent value
 */
function updateSubPlotHeight(val) {
    document.getElementById("chart-height-output").value = val + "%";
    var gridOption = chart.getOption().grid;
    var len = (amplicon || geneFeature) ? gridOption.length - 1 : gridOption.length;
    for (var i = 0; i < len; i++) { // Do not adjust gene feature height
        gridOption[i]["height"] = val + "%";
        if (i > 0) {
            // After adjusting height, need to adjust top margin as well
            gridOption[i]["top"] =
                parseFloat(gridOption[i - 1]["top"]) +
                parseFloat(gridOption[i - 1]["height"]) +
                parseFloat(document.getElementById("chart-top-input").value) + "%"; // set accoring to user's settings
        }
    };
    if (amplicon || geneFeature){
        gridOption[len]["top"] =
                parseFloat(gridOption[len - 1]["top"]) +
                parseFloat(gridOption[len - 1]["height"]) +
                parseFloat(document.getElementById("chart-top-input").value) + "%";
    }
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
                (parseFloat(gridOption[i - 1]["top"].replace("%","")) +
                parseFloat(gridOption[i - 1]["height"]) +
                parseFloat(val)).toFixed(1) + "%";
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
    var len = (amplicon || geneFeature) ? yAxisOption.length - 1 : yAxisOption.length
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
 * Reset Grid Dislay to optimal configuration
 */
function resetGridDisplay(){
    var chartOption = chart.getOption();
    var currentSamples = getCurrentSamples(chartOption);
    var gridOption = wgscovplot.getGrids(currentSamples, geneFeature, amplicon);
    chart.setOption({grid: gridOption});
    updateControlMenu();
}

/**
 * Dispatch click/dbclick actions for the whole chart
 */
function onChartDataZoomActions(){

    chart.on("click", function(params){
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === 'custom'){
            document.getElementById("start-pos").value = params.value.start;
            document.getElementById("end-pos").value = params.value.end;
            setDataZoom(params.value.start, params.value.end);
        }

    });

    chart.on("dblclick", function(params){
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === 'custom'){
            document.getElementById("start-pos").value = 1;
            document.getElementById("end-pos").value = refSeqLength;
            setDataZoom(1, refSeqLength);
        }
    });

}

/**
 * The Control Menu is updated when the number of selected samples changes
 * Menu is updated to reflect chart properties such as subplot height/top/left/right margin
 */
function updateControlMenu() {
    var gridOption = chart.getOption().grid;
    if (gridOption.length > 0) {
        var height = parseFloat(gridOption[0].height);
        var top = parseFloat(gridOption[0].top);
        var left = parseFloat(gridOption[0].left);
        var right = parseFloat(gridOption[0].right);
        document.getElementById("chart-height-input").value = height;
        document.getElementById("chart-height-output").value = height + "%";
        // Because height of each subplot = (1/(n+verticalRatio)) * 100 - heightOffset(6) + "%",
        // top margin is set 4 so need to plus 2.
        document.getElementById("chart-top-input").value = top + 2.0;
        document.getElementById("chart-top-output").value = top + 2.0 + "%";
        //update left
        document.getElementById("chart-left-input").value = left;
        document.getElementById("chart-left-output").value = left+ "%";

        document.getElementById("chart-right-input").value = right;
        document.getElementById("chart-right-output").value = right + "%";
        if (amplicon || geneFeature){
            document.getElementById("genefeature-height-input").value = parseFloat(gridOption[gridOption.length-1].height);
            document.getElementById("genefeature-height-output").value = parseFloat(gridOption[gridOption.length-1].height)+ "%";
        }
    }
}

