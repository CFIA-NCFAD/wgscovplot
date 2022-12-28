// Import the echarts core module, which provides the necessary interfaces for using echarts.
import * as echarts from "echarts/core";

import {BarChart, CustomChart, HeatmapChart, LineChart} from "echarts/charts";
import {
    DatasetComponent,
    DataZoomComponent,
    GridComponent,
    MarkAreaComponent,
    MarkLineComponent,
    TitleComponent,
    ToolboxComponent,
    TooltipComponent,
    VisualMapComponent
} from "echarts/components";
import {CanvasRenderer, SVGRenderer} from "echarts/renderers";
import {getCoverageChartOption, getGrids, getTooltips} from "./chart";

import {getCoverageThresholdLine, getMarkArea, getRegionAmpliconDepthRenderer} from "./chartOptions/series";
import {isNil, maxBy} from "lodash";

import $ from "jquery";
import {getGeneFeatureRenderer} from "./features";
import {getVariantHeatmapOption} from "./variantHeatmap";

echarts.use(
    [TooltipComponent, GridComponent,
        LineChart, BarChart, ToolboxComponent, HeatmapChart,
        DataZoomComponent, CustomChart, VisualMapComponent, TitleComponent,
        DatasetComponent, SVGRenderer, CanvasRenderer, MarkAreaComponent, MarkLineComponent]
);

/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are reserved
 * Users's settings are respected by keeping old settings and set it back.
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @param {Elements} elements - Object of HTML element names to HTMLElement object
 */
function updateCoverageChartOption({db, elements}) {
    const {
        chart,
        show_genes,
    } = db;
    let chartOption = chart.getOption();
    db.tooltipEnabled = elements.$toggleTooltip.checked;
    db.tooltipTriggerOn = (db.tooltipEnabled) ? (elements.$toggleShowTooltipOnClick.checked ? "click" : "mousemove") : "none";
    db.showVariantSiteTooltips = elements.$toggleShowVariantSiteTooltips.checked;
    db.showNonVariantSiteTooltips = elements.$toggleShowNonVariantSiteTooltips.checked;
    db.crossSampleComparisonInTooltips = elements.$toggleTooltipVariantComparison.checked;
    db.showCovStatsInTooltips = elements.$toggleTooltipCovStats.checked;
    db.showVariantLabels = elements.$toggleShowVariantLabels.checked;
    db.showGeneLabels = elements.$toggleXAxisLabel.checked;
    db.hideOverlappingVariantLabels = elements.$toggleHideOverlappingVariantLabels.checked;
    db.low_coverage_threshold = parseInt(elements.$lowThreshold.value);
    let updateOption = getCoverageChartOption(db);

    let lowCoverageRegion = elements.$toggleShowLowCovRegions.checked;
    db.showLowCovRegionsOpacity = lowCoverageRegion ? 0.4 : 0.0;

    // Preserve tooltip config in series options
    let seriesOption = updateOption.series;
    seriesOption.forEach(element => {
        if (element.type === "line") {
            element.markArea.itemStyle.opacity = db.showLowCovRegionsOpacity;
            element.tooltip.trigger = db.showNonVariantSiteTooltips ? "axis" : "none";
        } else if (element.type === "bar") {
            element.tooltip.trigger = db.showVariantSiteTooltips ? "axis" : "none";
        }
    });

    let fixTooltipPosition = document.getElementById("fix-tooltip-position").checked;
    updateOption.tooltip[0].position = tooltipPosition(fixTooltipPosition);

    // preserve DataZoom config
    updateOption.dataZoom = chartOption.dataZoom;
    updateOption.dataZoom.forEach(element => {
        element.xAxisIndex = [...Array(updateOption.grid.length).keys()];
    });

    // preserve y-axis limit/max value and scale type
    db.scaleType = elements.$selectYScale.value;
    db.yMax = parseInt(elements.$ymax.value);
    updateOption.yAxis = updateYAxisOption({yAxisOption: updateOption.yAxis, db});

    if (!isNil(db.amplicon_depths)) {
        db.show_amplicons = elements.$toggleShowAmplicon.checked;
        for (let i = 0; i < seriesOption.length - 1; i++) {
            if (seriesOption[i].type === "custom") {
                seriesOption[i].renderItem = getRegionAmpliconDepthRenderer(db.show_amplicons);
            }
        }
        updateOption.grid = getGrids(db);
    }

    if (show_genes) {
        seriesOption[seriesOption.length - 1].renderItem = getGeneFeatureRenderer(db);
    }

    db.leftMargin = elements.$chartLeftInput.value + "%";
    db.rightMargin = elements.$chartRightInput.value + "%";
    // Reserve grid option
    updateOption.grid.forEach(element => {
        element.left = db.leftMargin;
        element.right = db.rightMargin;
    });

    //set chart option
    chart.setOption(updateOption, {notMerge: true});
    // Update control menu
    updateControlMenu({db, elements});
    db.variantHeatmap.setOption(getVariantHeatmapOption(db));
}

/**
 * Set scale for y Axis
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @param {Elements} elements - Object of HTML element names to HTMLElement object
 */
function setScale({db, elements}) {
    db.scaleType = elements.$selectYScale.value;
    db.yMax = parseInt(elements.$ymax.value);
    let yAxisOption = db.chart.getOption().yAxis;
    let updatedYAxisOption = updateYAxisOption({yAxisOption, db})
    db.chart.setOption({yAxis: updatedYAxisOption});
}

/**
 * Set yMax for Y Axis
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @param {Elements} elements - Object of HTML element names to HTMLElement object
 */
function setYMax({db, elements}) {
    const {
        show_amplicons,
        show_genes,
    } = db;
    let yMax = elements.$ymax.value;
    let yAxisOption = db.chart.getOption().yAxis;
    let len = (show_amplicons || show_genes) ? yAxisOption.length - 1 : yAxisOption.length;
    yAxisOption.forEach(element => {
        if (element.gridIndex < len) {
            element.max = yMax;
        }
    });
    db.chart.setOption({yAxis: yAxisOption});
}

/**
 * Update scale type and max for Y Axis
 * @param {Object} yAxisOption - Options of Yaxis need to be updated
 * @param {Object} db - wgscovplot DB object
 * @returns {Object} Returns the updated options (Scale type or ymax) for yAxis
 */
function updateYAxisOption({yAxisOption, db}) {
    const {
        show_amplicons,
        show_genes,
        scaleType,
        yMax
    } = db;
    let len = (show_amplicons || show_genes) ? yAxisOption.length - 1 : yAxisOption.length;
    if (scaleType === "value") {
        // Linear scale
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 0;
                element.max = yMax;
            }
        });
    } else {
        // Log scale
        yAxisOption.forEach(element => {
            if (element.gridIndex < len) {
                element.type = scaleType;
                element.min = 1;
                element.max = yMax;
            }
        });
    }
    return yAxisOption;
}

/**
 * Get selected samples from sample select2 input
 *
 * If the chart is not initalizaed yet, get 3 first samples from list of all samples that could be displayed.
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 * @returns {string[]} Selected sample names
 */
function getSelectedSamples({db, elements}) {
    let selectedSamples = [];
    const chartOptions = db.chart.getOption();
    if (isNil(chartOptions)) {
        selectedSamples = db.samples.slice(0, 3);
    } else {
        let selectData = elements.$selectedSamples.select2("data");
        for (let {text} of Object.values(selectData)) {
            selectedSamples.push(text);
        }
    }
    return selectedSamples;
}

/**
 * Init chart display option event handlers for adjusting plot margins, heights, etc
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function initChartDisplayEventHandlers({db, elements}) {
    const {
        $chartTopInput,
        $chartRightInput,
        $chartHeightInput,
        $chartLeftInput,
        $renderEnv,
        $toggleDarkMode,
        $startPos,
        $endPos,
    } = elements;
    $startPos.value = 1;
    $endPos.value = db.positions.length;
    $chartTopInput.addEventListener(
        "change",
        () => updateSubPlotTopMargin({
            val: $chartTopInput.value,
            db,
            elements,
        })
    );
    $chartHeightInput.addEventListener(
        "change",
        () => updateSubPlotHeight({
            val: $chartHeightInput.value,
            db,
            elements,
        })
    );

    $chartLeftInput.addEventListener("change", () => {
        db.leftMargin = $chartLeftInput.value;
        elements.$chartLeftOutput.value = db.leftMargin + "%";
        let gridOption = db.chart.getOption().grid;
        gridOption.forEach(element => {
            element.left = db.leftMargin + "%";
        });
        db.chart.setOption({grid: gridOption});
    });

    $chartRightInput.addEventListener("change", () => {
        db.rightMargin = $chartRightInput.value
        elements.$chartRightOutput.value = db.rightMargin + "%";
        let gridOption = db.chart.getOption().grid;
        gridOption.forEach(element => {
            element.right = db.rightMargin + "%";
        });
        db.chart.setOption({grid: gridOption});
    });

    $renderEnv.addEventListener(
        "change",
        () => initWgscovplotRenderEnv({db, elements})
    );

    $toggleDarkMode.addEventListener(
        "change",
        () => initWgscovplotRenderEnv({db, elements})
    );
    if (db.show_genes === true) {
        elements.$genefeatureHeightInput.addEventListener("change", () => {
            const val = elements.$genefeatureHeightInput.value;
            elements.$genefeatureHeightOutput.value = val + "%";
            let gridOption = db.chart.getOption().grid;
            gridOption[gridOption.length - 1].height = val + "%";
            db.chart.setOption({grid: gridOption});
        });
    }


}

/**
 * Init chart tooltip event handlers
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function initTooltipEventHandlers({db, elements}) {
    const {
        chart
    } = db;

    elements.$toggleTooltip.addEventListener("change", () => {
        db.tooltipEnabled = elements.$toggleTooltip.checked;
        db.tooltipTriggerOn = db.tooltipEnabled ? (elements.$toggleShowTooltipOnClick.checked ? "click" : "mousemove") : "none";
        chart.setOption({
            tooltip: {triggerOn: db.tooltipTriggerOn},
        });
    });

    elements.$toggleFixedTooltipPosition.addEventListener("change", () => {
        db.fixedTooltipPosition = elements.$toggleFixedTooltipPosition.checked;
        let tooltip = chart.getOption().tooltip;
        tooltip[0].position = tooltipPosition(db.fixedTooltipPosition);
        chart.setOption({tooltip});
    });

    elements.$toggleShowTooltipOnClick.addEventListener("change", () => {
        db.tooltipEnabled = elements.$toggleTooltip.checked;
        db.tooltipTriggerOn = db.tooltipEnabled ? (elements.$toggleShowTooltipOnClick.checked ? "click" : "mousemove") : "none";
        chart.setOption({
            tooltip: {triggerOn: db.tooltipTriggerOn}
        });
    });

    elements.$toggleShowVariantSiteTooltips.addEventListener("change", () => {
        db.showVariantSiteTooltips = elements.$toggleShowVariantSiteTooltips.checked;
        updateTooltipOption({db, elements});
    });

    elements.$toggleShowNonVariantSiteTooltips.addEventListener("change", () => {
        db.showNonVariantSiteTooltips = elements.$toggleShowNonVariantSiteTooltips.checked;
        updateTooltipOption({db, elements});
    });

    elements.$toggleTooltipVariantComparison.addEventListener("change", () => {
        db.crossSampleComparisonInTooltips = elements.$toggleTooltipVariantComparison.checked;
        updateTooltipOption({db, elements});
    });

    elements.$toggleTooltipCovStats.addEventListener("change", () => {
        db.showCovStatsInTooltips = elements.$toggleTooltipCovStats.checked;
        updateTooltipOption({db, elements});
    });
}

/**
 *
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function initLowCovRegionHighlighting({db, elements}) {
    elements.$btnSetLowCovThreshold.addEventListener("click", () => {
        db.low_coverage_threshold = parseInt(elements.$lowThreshold.value);
        db.showLowCovRegions = elements.$toggleShowLowCovRegions.checked;
        db.showLowCovRegionsOpacity = db.showLowCovRegions ? 0.5 : 0.0;
        db.selectedSamples = getSelectedSamples({db, elements});
        let series = db.chart.getOption().series;
        for (let i = 0; i < db.selectedSamples.length; i++) {
            let sample = db.selectedSamples[i];
            series[i].markArea = getMarkArea({sample, db});
            series[i].markLine = getCoverageThresholdLine(db);
        }
        db.chart.setOption({series});
    })
}

/**
 * Initialize events handlers for ECharts and form control elements
 *
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @param {Elements} elements - jQuery selection objects or HTML control elements
 */
function initEventHandlers({db, elements}) {
    const {
        chart,
        ref_seq,
    } = db;
    const {
        $toggleAmplicons,
        $toggleShowVariantLabels,
        $toggleHideOverlappingVariantLabels,
        $chartLeftInput,
        $chartRightInput,
        $toggleXAxisLabel,
    } = elements;
    initChartDisplayEventHandlers({db, elements});
    initTooltipEventHandlers({db, elements});
    initLowCovRegionHighlighting({db, elements});

    /**
     * Events for selecting sample
     */
    elements.$selectedSamples.on("select2:select", (evt) => {
        let element = evt.params.data.element;
        let $element = $(element);
        $element.detach();
        elements.$selectedSamples.append($element);
        elements.$selectedSamples.trigger("change");
    });

    elements.$selectedSamples.on("change", () => {
        const selectedSamples = [];
        let selectData = elements.$selectedSamples.select2("data");
        for (let {text} of Object.values(selectData)) {
            selectedSamples.push(text);
        }
        db.selectedSamples = selectedSamples;
        updateCoverageChartOption({db, elements});
    });

    elements.$selectedGeneFeatures.select2({
        tags: true,
        width: "100%",
    });

    /**
     * Events for Axis Options
     */
    elements.$selectYScale.addEventListener(
        "change",
        () => {
            setScale({db, elements})
        })

    elements.$btnSetYMax.addEventListener(
        "click",
        () => {
            setYMax({db, elements})
        })

    /**
     * Display Options
     */
    elements.$btnShowZoomRegion.addEventListener(
        "click",
        () => {
            setDataZoom({
                start: 0,
                end: 0,
                db,
                elements,
            });
        }
    )

    elements.$btnResetView.addEventListener("click", () => {
        setDataZoom({
            start: 1,
            end: db.positions.length,
            db,
            elements
        })
    });

    if (db.show_genes === true) {
        elements.$toggleGeneLabel.addEventListener("change", () => {
            let seriesOption = db.chart.getOption().series;
            db.showGeneLabels = elements.$toggleGeneLabel.checked;
            seriesOption[seriesOption.length - 1].renderItem = getGeneFeatureRenderer(db);
            db.chart.setOption({series: [...seriesOption]});
        });
    }

    if (db.show_amplicons === true) {
        $toggleAmplicons.addEventListener("change", () => {

            db.show_amplicons = $toggleAmplicons.checked;
            let seriesOption = db.chart.getOption().series;
            for (let i = 0; i < seriesOption.length - 1; i++) {
                if (seriesOption[i].type === "custom") {
                    seriesOption[i].renderItem = getRegionAmpliconDepthRenderer(db.show_amplicons);
                }
            }
            db.showGeneLabels = elements.$toggleGeneLabel.checked;
            if (db.showGeneLabels) {
                seriesOption[seriesOption.length - 1].renderItem = getGeneFeatureRenderer(db);
            }
            const gridOptions = getGrids(db);
            db.leftMargin = $chartLeftInput.value;
            db.rightMargin = $chartRightInput.value;
            gridOptions.forEach(element => {
                element.left = db.leftMargin + "%";
                element.right = db.rightMargin + "%";
            });
            db.chart.setOption({series: [...seriesOption], grid: [...gridOptions]});
            updateControlMenu({db, elements});
        });
    }


    /*
    elements.$btnShowSelectedFeatures.addEventListener("click", () => {
        applyFeatureView({db, elements})
    });
    elements.$selectedGeneFeatures.on("change", function () {
        elements.$selectedGeneFeatures.select2("data");
    });
     */

    $toggleShowVariantLabels.addEventListener("change", () => {
        let series = db.chart.getOption().series;
        console.log(series)
        db.showVariantLabels = $toggleShowVariantLabels.checked;
        series.forEach(s => {
            if (s.type === "bar") {
                s.label.show = db.showVariantLabels;
            }
        });
        db.chart.setOption({series});
    });

    $toggleHideOverlappingVariantLabels.addEventListener("change", () => {
        let series = db.chart.getOption().series;
        db.hideOverlappingVariantLabels = $toggleHideOverlappingVariantLabels.checked;
        series.forEach(s => {
            if (s.type === "bar") {
                s.labelLayout.hideOverlap = db.hideOverlappingVariantLabels;
            }
        });
        db.chart.setOption({series});
    });

    $toggleXAxisLabel.addEventListener("change", () => {
        let xAxis = db.chart.getOption().xAxis;
        let showAxisLabel = $toggleXAxisLabel.checked;
        let gridLength = (db.show_amplicons || db.show_genes) ? xAxis.length - 1 : xAxis.length;
        for (let i = 0; i < gridLength; i++) {
            xAxis[i].axisLabel.show = showAxisLabel;
        }
        db.chart.setOption({xAxis});
    });

    // elements.$toggleShowLowCovRegionLabels.addEventListener("change", () => {
    //     let series = chart.getOption().series;
    //     let showLabel = elements.$toggleShowLowCovRegionLabels.checked;
    //     series.forEach(series => {
    //         if (series.type === "line") {
    //             series.markArea.label.show = showLabel;
    //         }
    //     });
    //     chart.setOption({series});
    // });

    elements.$toggleShowLowCovRegions.addEventListener("change", () => {
        let series = db.chart.getOption().series;
        db.showLowCovRegions = elements.$toggleShowLowCovRegions.checked;
        db.showLowCovRegionsOpacity = db.showLowCovRegions ? 0.4 : 0;
        elements.$btnSetLowCovThreshold.disabled = !db.showLowCovRegions;
        series.forEach(s => {
            if (s.type === "line") {
                s.markArea.itemStyle.opacity = db.showLowCovRegionsOpacity;
                s.markLine.lineStyle.opacity = db.showLowCovRegionsOpacity;
            }
        });
        db.chart.setOption({series});
    });

    /**
     * Toggle slider zoom
     */

    elements.$toggleDataZoomSlider.addEventListener("change", () => {
        db.showDataZoomSlider = elements.$toggleDataZoomSlider.checked;
        let nGrids = chart.getOption().grid.length;
        let xAxisIndex = [...Array(nGrids).keys()];
        chart.setOption({
            dataZoom: [
                {
                    type: "inside",
                    filterMode: "none",
                    xAxisIndex,
                },
                {
                    type: "slider",
                    show: db.showDataZoomSlider,
                    filterMode: "none",
                    xAxisIndex: db.showDataZoomSlider ? xAxisIndex : null,
                },
            ],
        });
        /*
        db.variantHeatmap.setOption({
            dataZoom: [
                {
                    type: "inside"
                },
                {
                    type: "slider",
                    show: db.showDataZoomSlider
                }
            ]
        });*/
    });

    elements.$btnResetMargins.addEventListener("click", () => {
        resetGridDisplay({db, elements});
    })

    /*
    elements.$btnSetYMax.addEventListener("click", () => {
        db.yMax = parseInt(elements.$ymax.value);
        let yAxisOption = db.chart.getOption().yAxis;
        yAxisOption = updateYAxisOption({yAxisOption, db});
        db.chart.setOption({yAxis: yAxisOption});
    });

    elements.$selectYScale.addEventListener("change", () => {
        db.yMax = maxBy(Object.values(db.mosdepth_info), "max_depth").max_depth;
        elements.$ymax.value = db.yMax;
        db.scaleType = elements.$selectYScale.value;
        let yAxisOption = db.chart.getOption().yAxis;
        yAxisOption = updateYAxisOption({yAxisOption, db});
        db.chart.setOption({yAxis: yAxisOption});
    });
     */

}

/**
 * Initialize the environment for the chart (sgv/canvas or dark/white mode)
 * The entire chart will be disposed and re-initialized
 * However, the old settings of charts are reserved (users's settings are respected)
 */
function initWgscovplotRenderEnv({db, elements}) {
    const {
        chart,
        variantHeatmap,
    } = db;
    const {
        $chart,
        $renderEnv,
        $toggleDarkMode,
    } = elements;
    let chartOptions = chart.getOption();
    if (isNil(chartOptions)) {
        chart.setOption(getCoverageChartOption(db));
        //variantHeatmap.setOption(getVariantHeatmapOption(db));
    } else {
        let renderEnv = $renderEnv.value;
        let isChecked = $toggleDarkMode.checked;
        let mode = isChecked ? "dark" : "white";
        let gridOption = chartOptions.grid;
        let seriesOption = chartOptions.series;
        let dataZoomOption = chartOptions.dataZoom;
        let tooltipOption = chartOptions.tooltip;
        // chart.dispose()
        echarts.dispose(chart); // destroy chart instance and re-init chart
        db.chart = echarts.init($chart, mode, {renderer: renderEnv});
        let option = getCoverageChartOption(db);
        // Keep grid option
        option.grid = gridOption;
        // Keep data zoom option
        option.dataZoom = dataZoomOption;
        // Keep yAxis option
        option.yAxis = updateYAxisOption({yAxisOption: option.yAxis, db});
        // Keep tooltip
        option.tooltip = tooltipOption;
        // Keep series
        option.series = seriesOption;
        //set chart option
        db.chart.setOption(option);
    }
    updateControlMenu({db, elements});
    onChartDataZoomActions({db, elements});
}

/**
 * Update option for displaying tooltip for variant and non variant sites
 * @param {Object} db - An array of samples name
 * @param {Object} elements -
 */
function updateTooltipOption(
    {
        db,
        elements,
    }) {
    const {
        chart,
        showVariantSiteTooltips,
        showNonVariantSiteTooltips,
    } = db;
    db.tooltipEnabled = elements.$toggleTooltip.checked;
    db.tooltipTriggerOn = db.tooltipEnabled ? (elements.$toggleShowTooltipOnClick.checked ? "click" : "mousemove") : "none";
    let series = chart.getOption().series;
    series.forEach(element => {
        if (element.type === "line") {
            element.tooltip.trigger = showNonVariantSiteTooltips ? "axis" : "none";
        } else if (element.type === "bar") {
            element.tooltip.trigger = showVariantSiteTooltips ? "axis" : "none";
        }
    });
    db.showCovStatsInTooltips = elements.$toggleTooltipCovStats.checked;
    let tooltip = getTooltips(db);
    db.fixedTooltipPosition = elements.$toggleFixedTooltipPosition.checked;
    tooltip[0].position = tooltipPosition(db.fixedTooltipPosition);
    tooltip[0].triggerOn = db.tooltipTriggerOn;
    chart.setOption({tooltip, series});
}

/**
 * Functions to return tooltip position
 * @param {boolean} isChecked
 * @returns Returns the tooltip position
 */
function tooltipPosition(isChecked) {
    if (isChecked) {
        return function (pos, params, dom, rect, size) {
            let obj = {
                top: 5,
            };
            obj[(pos[0] < size.viewSize[0] / 2) ? "right" : "left"] = 5;
            return obj;
        };
    } else {
        return "cursor"; // follow cursor
    }
}

/**
 * Adjust Variant Heatmap height
 * @param {number} val - Subplots height percent value
 */
function updateVarMapHeight(val) {
    document.getElementById("varmap-height-output").value = val;
    let gridOption = variantHeatmap.getOption().grid;
    gridOption[0].height = val + "%";
    variantHeatmap.setOption({grid: gridOption});
}

/**
 * Apply View for selected gene feature
 */
function applyFeatureView({db, elements}) {
    const {
        echart_features,
        ref_seq,
    } = db;
    const {
        $select2Features,
    } = elements;
    let features = $select2Features.select2("data");
    let starts = [], ends = [];
    let minStart, maxEnd;
    features.forEach(({text: feature}) => {
        for (let {name, value: {start, end}} of echart_features) {
            if (feature === name) {
                starts.push(start);
                ends.push(end);
                break;
            }
        }
    });
    minStart = Math.min(...starts);
    maxEnd = Math.max(...ends);
    if (features.length === 0) {
        minStart = 1;
        maxEnd = ref_seq.length;
    }
    setDataZoom({
        start: minStart,
        end: maxEnd,
        db,
        elements,
    });
}

/**
 * Adjust subplot heights and top margins and chart-height-output value.
 *
 * Triggered on change of #chart-height-input value.
 *
 * @param {string} val - Subplots height percent value
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function updateSubPlotHeight(
    {
        val,
        db,
        elements,
    }) {
    const {
        chart,
        show_amplicons,
        show_genes,
    } = db;
    const {
        $chartHeightOutput,
        $chartTopInput,
    } = elements;
    $chartHeightOutput.value = val + "%";
    let gridOption = chart.getOption().grid;
    let len = (show_amplicons || show_genes) ? gridOption.length - 1 : gridOption.length;
    let topInputValue = parseFloat($chartTopInput.value);
    for (let i = 0; i < len; i++) {
        let grid = gridOption[i];
        grid.height = val + "%";
        if (i > 0) {
            // After adjusting height, need to adjust top margin as well
            let gridPrevious = gridOption[i - 1];
            let top = parseFloat(gridPrevious.top);
            let height = parseFloat(gridPrevious.height);
            grid.top = top + height + topInputValue + "%";
        }
    }
    if (show_amplicons || show_genes) {
        let gridLastPlot = gridOption[len - 1];
        let top = parseFloat(gridLastPlot.top);
        let height = parseFloat(gridLastPlot.height);
        gridOption[len].top = top + height + topInputValue + "%";
    }
    chart.setOption({grid: gridOption});
}


/**
 * Adjust top margin of chart
 * @param {WgsCovPlotDB} db
 * @param {string} val - Subplots top margin percent value
 * @param {HTMLOutputElement} $chartTopOutput - chart right margin output HTML element
 */
function updateSubPlotTopMargin(
    {
        val,
        db,
        elements: {
            $chartTopOutput,
        }
    }) {
    const {
        chart,
    } = db;
    $chartTopOutput.value = val + "%";
    let gridOption = chart.getOption().grid;
    for (let i = 0; i < gridOption.length; i++) {
        let grid = gridOption[i];
        if (i === 0) {
            grid.top = val + "%";
        } else {
            let previousGrid = gridOption[i - 1];
            let previousTop = parseFloat(previousGrid.top.replace("%", ""));
            let previousHeight = parseFloat(previousGrid.height);
            let topValue = parseFloat(val);
            grid.top = (previousTop + previousHeight + topValue).toFixed(1) + "%";
        }
    }
    chart.setOption({grid: gridOption});
}

/**
 * Set zoom view for the chart
 * @param {number} start - Start view point
 * @param {number} end - End view point
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function setDataZoom(
    {
        start = 0,
        end = 0,
        db,
        elements,
    }) {
    const {
        $startPos,
        $endPos,
    } = elements;
    if (start === 0 && end === 0) {
        start = parseInt($startPos.value);
        end = parseInt($endPos.value);
    } else {
        // update view range input element values
        $startPos.value = `${start}`;
        $endPos.value = `${end}`;
    }
    db.chart.dispatchAction({
        type: "dataZoom",
        startValue: start,
        endValue: end,
    });
}

/**
 * Reset Grid Dislay to optimal configuration
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function resetGridDisplay({db, elements}) {
    if (db.show_amplicons === true){
        db.show_amplicons = elements.$toggleAmplicons.checked;
    }
    let grid = getGrids(db);
    db.chart.setOption({grid});
    updateControlMenu({db, elements});
}

/**
 * Dispatch click/dbclick actions for the whole chart
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function onChartDataZoomActions({db, elements}) {
    const {chart, ref_seq} = db;
    const {$startPos, $endPos} = elements;

    chart.on("click", function (
        {
            componentIndex,
            componentSubType,
            value: {start, end},
        }) {
        if (componentIndex === chart.getOption().series.length - 1 && componentSubType === "custom") {
            $startPos.value = start;
            $endPos.value = end;
            setDataZoom({
                start,
                end,
                db,
                elements,
            });
        }
    });

    chart.on("dblclick", function ({componentIndex, componentSubType}) {
        if (componentIndex === chart.getOption().series.length - 1 && componentSubType === "custom") {
            const start = 1;
            const end = ref_seq.length;
            $startPos.value = start;
            $endPos.value = end;
            setDataZoom({
                start,
                end,
                db,
                elements,
            });
        }
    });
}

/**
 * The Control Menu is updated when the number of selected samples changes
 * Menu is updated to reflect chart properties such as subplot height/top/left/right margin
 * @param {WgsCovPlotDB} db
 * @param {Elements} elements
 */
function updateControlMenu({db, elements}) {
    let gridOption = db.chart.getOption().grid;
    if (gridOption.length > 0) {
        let height = parseFloat(gridOption[0].height);
        let top = parseFloat(gridOption[0].top);
        let left = parseFloat(gridOption[0].left);
        let right = parseFloat(gridOption[0].right);
        elements.$chartHeightInput.value = height;
        elements.$chartHeightOutput.value = height + "%";
        // Because height of each subplot = (1/(n+verticalRatio)) * 100 - heightOffset(6) + "%",
        // top margin is set 4 so need to plus 2.
        elements.$chartTopInput.value = top + 2.0;
        elements.$chartTopOutput.value = top + 2.0 + "%";
        // update left margin
        elements.$chartLeftInput.value = left;
        elements.$chartLeftOutput.value = left + "%";
        // update right margin
        elements.$chartRightInput.value = right;
        elements.$chartRightOutput.value = right + "%";
        if (db.show_amplicons || db.show_genes) {
            let featuresGrid = gridOption[gridOption.length - 1];
            elements.$genefeatureHeightInput.value = featuresGrid.height;
            elements.$genefeatureHeightOutput.value = featuresGrid.height;
        }
    }
}

export {
    initEventHandlers,
    initWgscovplotRenderEnv,
    echarts,
}
