/**
 * Initialize the events handler for chart
 */
function initWgscovplotEvent() {
    $(document).ready(function () {

        /**
         * Toggle dark mode, the entire chart will be disposed and re-initialized
         * However, the old settings of charts are set back (users's settings are respected)
         */
        $("#toggle-darkmode").change(function () {
            initWgscovplotRenderEnv();
        });

        $("#selected-samples").select2({
            tags: true,
        });

        /**
         * Jquery actions to make the list of samples is not forced in alphabetical order
         */
        $("#selected-samples").on("select2:select", function (evt) {
            let element = evt.params.data.element;
            let $element = $(element);
            $element.detach();
            $(this).append($element);
            $(this).trigger("change");
        });

        $("#selected-segments").select2({
            tags: true,
        });

        /**
         * Change samples for viewing
         */
        $("#selected-samples").on("change", function () {
            let [samples, segments] = getCurrentSamplesSegments(chart.getOption());
            if (samples.length > 0 && segments.length > 0) {
                updateFluCoverageChartOption(samples, segments);
            }
        });

        /**
         * Change segments for viewing
         */
        $("#selected-segments").on("change", function () {
            let [samples, segments] = getCurrentSamplesSegments(chart.getOption());
            if (samples.length > 0 && segments.length > 0) {
                updateFluCoverageChartOption(samples, segments);
            }
        });

        /**
         * Toggle to show X Axis Label
         */
        $("#toggle-xaxis-label").change(function () {
            let xAxisOption = chart.getOption().xAxis;
            let showAxisLabel = $(this).prop("checked");
            for (let i = 0; i < xAxisOption.length - 1; i++) {
                xAxisOption[i].axisLabel.show = showAxisLabel;
            }
            chart.setOption({xAxis: [...xAxisOption]});
        });

        /**
         * Toggle tooltip for coverage chart
         */
        $("#toggle-tooltip").change(function () {
            let isChecked = $(this).prop("checked");
            let triggerType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
            chart.setOption({
                tooltip: {triggerOn: isChecked ? triggerType : "none"},
            });
        });

        /**
         * Fix tooltip position on the charts
         * If it is fixed: tooltip will be on the right if mouse hovering on the left and vice versa
         * Default tooltip follows cursor
         */
        $("#fix-tooltip-position").change(function () {
            let isChecked = $(this).prop("checked");
            let tooltipOption = chart.getOption().tooltip;
            tooltipOption[0].position = tooltipPosition(isChecked);
            chart.setOption({tooltip: tooltipOption});
        });

        /**
         * Toogle tooltip trigger on click
         */
        $("#toggle-tooltip-trigger-click").change(function () {
            let isChecked = $(this).prop("checked");
            let isTooltipEnable = document.getElementById("toggle-tooltip").checked;
            if (isTooltipEnable) {
                chart.setOption({
                    tooltip: {triggerOn: isChecked ? "click" : "mousemove"}
                });
            } else {
                chart.setOption({
                    tooltip: {triggerOn: "none"}
                });
            }
        });

        /**
         * Toggle tooltip display for variant sites
         */
        $("#toggle-tooltip-variant-sites").change(function () {
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let [samples, segments] = getCurrentSamplesSegments(chartOption);
            updateTooltipOption(samples, segments, window.depths, window.variants, seriesOption,
                isChecked, document.getElementById("toggle-tooltip-non-variant-sites").checked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });

        /**
         * Toggle tooltip display for non-variant sites
         */
        $("#toggle-tooltip-non-variant-sites").change(function () {
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let [samples, segments] = getCurrentSamplesSegments(chartOption);
            updateTooltipOption(samples, segments, window.depths, window.variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked, isChecked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });

        /**
         * Toggle tooltip display for variant comparison
         */
        $("#toggle-variant-comparison").change(function () {
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let [samples, segments] = getCurrentSamplesSegments(chartOption);
            updateTooltipOption(samples, segments, window.depths, window.variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked,
                document.getElementById("toggle-tooltip-non-variant-sites").checked,
                isChecked, document.getElementById("toggle-coverage-stat").checked);
        });

        /**
         * Toggle tooltip display for coverage statistics viewing
         */
        $("#toggle-coverage-stat").change(function () {
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let [samples, segments] = getCurrentSamplesSegments(chartOption);
            updateTooltipOption(samples, segments, window.depths, window.variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked,
                document.getElementById("toggle-tooltip-non-variant-sites").checked,
                document.getElementById("toggle-variant-comparison").checked, isChecked);
        });

        /**
         * Toggle to show Mutation below Variant Sites
         */
        $("#toggle-mutation").change(function () {
            let seriesOption = chart.getOption().series;
            let showMutation = $(this).prop("checked");
            seriesOption.forEach(series => {
                if (series.type === "bar") {
                    series.label.show = showMutation;
                }
            });
            chart.setOption({series: [...seriesOption]});
        });

        //toggle-low-coverage-regions-label
        $("#toggle-low-coverage-regions-label").change(function () {
            let seriesOption = chart.getOption().series;
            let showLabel = $(this).prop("checked");
            seriesOption.forEach(series => {
                if (series.type === "line") {
                    series.markArea.label.show = showLabel;
                }
            });
            chart.setOption({series: [...seriesOption]});
        });

        //toggle-low-coverage-regions
        $("#toggle-low-coverage-regions").change(function () {
            let seriesOption = chart.getOption().series;
            let isShow = $(this).prop("checked");
            let opacity;
            if (isShow === true){
                opacity = 0.4;
                document.getElementById("set-low-threshold").disabled = false;
            } else {
                opacity = 0.0;
                document.getElementById("set-low-threshold").disabled = true;
            }
            seriesOption.forEach(series => {
                if (series.type === "line") {
                    series.markArea.itemStyle.opacity = opacity;
                    series.markLine.lineStyle.opacity = opacity;
                }
            });
            chart.setOption({series: [...seriesOption]});
        });

        /**
         * Toggle to hide overlapping mutation under Variant Sites
         */
        $("#toggle-hideoverlap-mutation").change(function () {
            let seriesOption = chart.getOption().series;
            let isHideOverlap = $(this).prop("checked");
            seriesOption.forEach(series => {
                if (series.type === "bar") {
                    series.labelLayout.hideOverlap = isHideOverlap;
                }
            });
            chart.setOption({series: [...seriesOption]});
        });

    });
}

/**
 * Update option for displaying tooltip for variant and non-variant sites
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of segments
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The object of variants data
 * @param {Array<Array<Object>>} seriesOption - Series option of the chart
 * @param {boolean} variantSites - Whether to enable tooltip for variant sites or not
 * @param {boolean} nonVariantSites - Whether to enable tooltip for non-variant sites or not
 * @param {boolean} infoComparison - Whether to compare variant/ Coverage Stat information across selected samples
 * @param {boolean} coverageStatView - Whether to display Coverage Stat
 */
function updateTooltipOption(samples, segments, depths, variants, seriesOption,
                             variantSites, nonVariantSites, infoComparison, coverageStatView) {

    let tooltipEnable = document.getElementById("toggle-tooltip").checked;
    let triggerOnType;
    if (tooltipEnable) {
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    } else {
        triggerOnType = "none";
    }
    seriesOption.forEach(element => {
        if (element.type === "line") {
            element.tooltip.trigger = nonVariantSites ? "axis" : "none";
        } else if (element.type === "bar") {
            element.tooltip.trigger = variantSites ? "axis" : "none";
        }
    });
    let tooltipOption = wgscovplot.getFluTooltips(samples, segments, depths, variants, window.refSeq, window.refID,
        triggerOnType, infoComparison, coverageStatView, window.lowCoverageThreshold, window.primerData);
    let fixTooltipPosition = document.getElementById("fix-tooltip-position").checked;
    tooltipOption[0].position = tooltipPosition(fixTooltipPosition);
    tooltipOption[0].triggerOn = triggerOnType;
    chart.setOption({tooltip: tooltipOption, series: seriesOption});
}

/**
 * Functions to return tooltip position
 * @param {boolean} isChecked
 * @returns Returns the tooltip position
 */
function tooltipPosition(isChecked) {
    if (isChecked) {
        return function (pos, params, dom, rect, size) {
            let obj = {top: 5};
            obj[["left", "right"][+(pos[0] < size.viewSize[0] / 2)]] = 5;
            return obj;
        };
    } else {
        return "cursor"; // follow cursor
    }
}

/**
 * Update Coverage Chart option
 * @param {Array<string>} samples
 * @param {Array<string>} segments
 * @returns Returns the tooltip position
 */
function updateFluCoverageChartOption(samples, segments) {
    let scaleType;
    let chartOption = chart.getOption();
    let tooltipEnable = document.getElementById("toggle-tooltip").checked;
    let triggerOnType;
    if (tooltipEnable) {
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    } else {
        triggerOnType = "none";
    }
    let variantSites = document.getElementById("toggle-tooltip-variant-sites").checked;
    let nonVariantSites = document.getElementById("toggle-tooltip-non-variant-sites").checked;
    let variantComparison = document.getElementById("toggle-variant-comparison").checked;
    let coverageStatView = document.getElementById("toggle-coverage-stat").checked;
    let showMutation = document.getElementById("toggle-mutation").checked;
    let hideOverlapMutation = document.getElementById("toggle-hideoverlap-mutation").checked;
    let showXAxisLabel = document.getElementById("toggle-xaxis-label").checked;
    let threshold = document.getElementById("low-threshold").value;
    let updateOption = wgscovplot.getFluCoverageChartOption(samples, segments, window.depths, window.variants,
        window.refSeq, window.refID, window.lowCoverageRegions, threshold, window.primerData,
        triggerOnType, variantSites, nonVariantSites, variantComparison, coverageStatView, showMutation,
        showXAxisLabel, hideOverlapMutation);

    let lowCoverageRegion = document.getElementById("toggle-low-coverage-regions").checked;
    let opacity;
    if (lowCoverageRegion === true) {
        opacity = 0.4;
    } else {
        opacity = 0.0;
    }
    let seriesOption = updateOption.series;
    seriesOption.forEach(element => {
        if (element.type === "line") {
            element.markArea.itemStyle.opacity = opacity;
            element.tooltip.trigger = nonVariantSites ? "axis" : "none";
        } else if (element.type === "bar") {
            element.tooltip.trigger = variantSites ? "axis" : "none";
        }
    });

    let fixTooltipPosition = document.getElementById("fix-tooltip-position").checked;
    updateOption.tooltip[0].position = tooltipPosition(fixTooltipPosition);

    // Reserve datzoom
    let oldDataZoom = chartOption.dataZoom;
    updateOption.dataZoom = oldDataZoom;
    updateOption.dataZoom.forEach(element => {
        element.xAxisIndex = [...Array(updateOption.grid.length).keys()];
    });

    //Reserve scale and yAxis max
    scaleType = $("#scale").val();
    updateOption.yAxis = updateYAxisOption(updateOption.yAxis, scaleType);

    // Reserve grid option
    updateOption.grid.forEach(element => {
        element.left = $("#chart-left-input").val() + "%";
        element.right = $("#chart-right-input").val() + "%";
    });

    chart.setOption(option = updateOption, {notMerge: true});
    document.getElementById("ymax").value = chart.getOption().yAxis[0].max;
    updateControlMenu();
}

/**
 * Update scale type and max for Y Axis
 * @param {Object} yAxisOption - Options of Yaxis need to be updated
 * @param {string} scaleType - Either log or value
 * @returns {Object} Returns the updated options (Scale type or ymax) for yAxis
 */
function updateYAxisOption(yAxisOption, scaleType) {
    let len = yAxisOption.length - 1;
    if (scaleType === "value") {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 0;
            }
        });
    } else {
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 1;
            }
        });
    }
    return yAxisOption;
}

/**
 * Select default samples, segments when chart is initialized (the first 3 samples, segments are plotted)
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<string>} segments - An array of samples name
 */
function setDefaultSamplesSegments(samples, segments) {
    // Set default samples display
    let $selectedSamples = $("#selected-samples");
    $selectedSamples.select2();
    $selectedSamples.val(samples);
    $selectedSamples.trigger("change");

    // Set default samples display
    let $selectedSegments = $("#selected-segments");
    $selectedSegments.select2();
    $selectedSegments.val(segments);
    $selectedSegments.trigger("change");
}

/**
 * Get the list of current samples, segments on the control menu
 * @param {Object} chartOption
 * @return <Array<Array<string>>
 */
function getCurrentSamplesSegments(chartOption) {
    let samples = [];
    let segments = [];
    // If the chart is not initialized yet, get 3 first samples from window.samples
    if (chartOption === undefined || chartOption === null) {
        samples = window.samples.slice(0, 3);
        segments = window.segments.slice(0, 8);
    } else {
        let selectData1 = $("#selected-samples").select2("data");
        let selectData2 = $("#selected-segments").select2("data");
        for (let [key1, entries1] of selectData1.entries()) {
            samples.push(selectData1[key1].text);
        }
        for (let [key2, entries2] of selectData2.entries()) {
            segments.push(selectData2[key2].text);
        }
    }
    return [samples, segments];
}

/**
 * Initialize Chart Environment
 */
function initWgscovplotRenderEnv() {

    let chartOption = chart.getOption();
    const [plotSamples, plotSegments] = getCurrentSamplesSegments(chartOption);
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamplesSegments(plotSamples, plotSegments);
        chart.setOption(option = wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, window.depths,
            window.variants, window.refSeq, window.refID, window.lowCoverageRegions, window.lowCoverageThreshold, window.primerData));
    } else {
        let renderEnv = document.getElementById("render-env").value;
        let isChecked = document.getElementById("toggle-darkmode").checked;
        let mode = isChecked ? "dark" : "white";
        let gridOption = chart.getOption().grid;
        let seriesOption = chart.getOption().series;
        let scaleType = chart.getOption().yAxis[0].type;
        //let yAxisMax = chart.getOption().yAxis[0].max;
        let dataZoomOption = chart.getOption().dataZoom;
        let tooltipOption = chart.getOption().tooltip;
        let lowCoverageRegion = document.getElementById("toggle-low-coverage-regions").checked;
        let opacity;
        if (lowCoverageRegion === true) {
            opacity = 0.4;
        } else {
            opacity = 0.0;
        }
        wgscovplot.echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        chart = wgscovplot.echarts.init($chart, mode, {renderer: renderEnv});
        let option = wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, window.depths, window.variants,
            window.refSeq, window.refID, window.lowCoverageRegions, window.lowCoverageThreshold, window.primerData);
        // Keep grid option
        option.grid = gridOption;
        // Keep data zoom option
        option.dataZoom = dataZoomOption;
        // Keep yAxis option
        option.yAxis = updateYAxisOption(option.yAxis, scaleType);
        // Keep tooltip
        option.tooltip = tooltipOption;
        // Keep series
        option.series = seriesOption;
        // Keep low coverage regions hightlight
        option.series.forEach(element => {
            if (element.type === "line") {
                element.markArea.itemStyle.opacity = opacity;
            }
        });
        //set chart option
        chart.setOption(option);
    }
    // Update yAxisMax on Control Menu:
    document.getElementById("ymax").value = chart.getOption().yAxis[0].max;
    onChartDataZoomActions();
}

/**
 * Set scale for y Axis
 */
function setScale() {
    let scaleType = document.getElementById("scale").value;
    let yAxisOption = updateYAxisOption(chart.getOption().yAxis, scaleType);
    chart.setOption({yAxis: yAxisOption});
}

/**
 * Set yMax for Y Axis
 */
function setYMax() {
    let yMax = document.getElementById("ymax").value;
    let yAxisOption = chart.getOption().yAxis;
    yAxisOption.forEach(element => {
        if (element.gridIndex < yAxisOption.length - 1) {
            element.max = yMax;
        }
    });
    chart.setOption({yAxis: yAxisOption});
}

/**
 * Set low coverage threshold with horizontal line in the charts
 */
function setLowCoverageThreshold() {
    let threshold = document.getElementById("low-threshold").value;
    let chartOption = chart.getOption();
    const [currentSamples, currentSegments] = getCurrentSamplesSegments(chartOption);
    let maxSegmentsLength = wgscovplot.getMaxSegmentsLength(currentSamples, currentSegments, window.depths);
    let segmentsInterval = wgscovplot.getSegmentsInterval(maxSegmentsLength);
    let seriesOption = chartOption.series;
    let lowCoverageRegion = document.getElementById("toggle-low-coverage-regions").checked;
    let opacity;
    if (lowCoverageRegion === true) {
        opacity = 0.4;
    } else {
        opacity = 0.0;
    }
    for (let i = 0; i < currentSamples.length; i++) {
        seriesOption[i].markArea = wgscovplot.getMarkArea(currentSamples[i], currentSegments, segmentsInterval, window.depths, parseInt(threshold), opacity);
        seriesOption[i].markLine = wgscovplot.getCoverageThresholdLine(segmentsInterval, threshold);
    }
    chart.setOption({series: seriesOption});
}

/**
 * Dispatch click/dbclick actions for the whole chart
 */
function onChartDataZoomActions() {

    chart.on("click", function (params) {
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === "custom") {
            setDataZoom(params.value.start, params.value.end);
        }

    });

    chart.on("dblclick", function (params) {
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === "custom") {
            let chartOption = chart.getOption();
            const [plotSamples, plotSegments] = getCurrentSamplesSegments(chartOption);
            const maxSegmentsLength = wgscovplot.getMaxSegmentsLength(plotSamples, plotSegments, window.depths);
            let xAxisMax = 1;
            xAxisMax = maxSegmentsLength.reduce((a, b) => {
                return a + b;
            });
            setDataZoom(1, xAxisMax);
        }
    });

}

/**
 * Set zoom view for the chart
 * @param {number} zoomStart - Start view point
 * @param {number} zoomEnd - End view point
 */
function setDataZoom(zoomStart, zoomEnd) {
    chart.dispatchAction({
        type: "dataZoom",
        startValue: zoomStart,
        endValue: zoomEnd,
    });
}

/**
 * Adjust subplot height
 * @param {number} val - Subplots height percent value
 */
function updateSubPlotHeight(val) {
    document.getElementById("chart-height-output").value = val + "%";
    let gridOption = chart.getOption().grid;
    let len = gridOption.length - 1;
    for (let i = 0; i < len; i++) { // Do not adjust gene feature height
        gridOption[i].height = val + "%";
        if (i > 0) {
            // After adjusting height, need to adjust top margin as well
            gridOption[i].top =
                parseFloat(gridOption[i - 1].top) +
                parseFloat(gridOption[i - 1].height) +
                parseFloat(document.getElementById("chart-top-input").value) + "%"; // set accoring to user's settings
        }
    }
    gridOption[len].top =
        parseFloat(gridOption[len - 1].top) +
        parseFloat(gridOption[len - 1].height) +
        parseFloat(document.getElementById("chart-top-input").value) + "%";

    chart.setOption({grid: gridOption});
}

/**
 * Adjust left margin of chart
 * @param {number} val - Left margin percent value
 */
function updateChartLeftMargin(val) {
    document.getElementById("chart-left-output").value = val + "%";
    let gridOption = chart.getOption().grid;
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
    let gridOption = chart.getOption().grid;
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
    let gridOption = chart.getOption().grid;
    for (let i = 0; i < gridOption.length; i++) {
        if (i === 0) {
            gridOption[i].top = val + "%";
        } else {
            gridOption[i].top =
                (parseFloat(gridOption[i - 1].top.replace("%", "")) +
                    parseFloat(gridOption[i - 1].height) +
                    parseFloat(val)).toFixed(1) + "%";
        }
    }
    chart.setOption({grid: gridOption});
}

/**
 * Reset Grid Display to optimal configuration
 */
function resetGridDisplay() {
    let chartOption = chart.getOption();
    const [currentSamples, currentSegments] = getCurrentSamplesSegments(chartOption);
    let gridOption = wgscovplot.getGrids(currentSamples, true, false, false);
    chart.setOption({grid: gridOption});
    updateControlMenu();
}

/**
 * The Control Menu is updated when the number of selected samples changes
 * Menu is updated to reflect chart properties such as subplot height/top/left/right margin
 */
function updateControlMenu() {
    let gridOption = chart.getOption().grid;
    if (gridOption.length > 0) {
        let height = parseFloat(gridOption[0].height);
        let top = parseFloat(gridOption[0].top);
        let left = parseFloat(gridOption[0].left);
        let right = parseFloat(gridOption[0].right);
        document.getElementById("chart-height-input").value = height;
        document.getElementById("chart-height-output").value = height + "%";
        // Because height of each subplot = (1/(n+verticalRatio)) * 100 - heightOffset(6) + "%",
        // top margin is set 4 so need to plus 2.
        document.getElementById("chart-top-input").value = top + 2.0;
        document.getElementById("chart-top-output").value = top + 2.0 + "%";
        //update left
        document.getElementById("chart-left-input").value = left;
        document.getElementById("chart-left-output").value = left + "%";

        document.getElementById("chart-right-input").value = right;
        document.getElementById("chart-right-output").value = right + "%";
    }
}
