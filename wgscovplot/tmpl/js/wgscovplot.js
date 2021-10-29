/**
 * Define the points of gene/amplicon features shape
 * @param {number} x - x-axis coordinate
 * @param {number} y  - x-axis coordinate
 * @param {number} width - width of shape
 * @param {number} height - height of shape
 * @param {number} strand - strand of feature
 * @param {string} feature - gene feature or amplicon feature
 * @returns {Array<Array<number>>} - Coordinate of 5 points for gene feature or 4 points for amplicon feature
 */
function renderPoints(x, y, width, height, strand, feature) {
    if (feature === 'gene_feature') {
        if (strand === 1) {
            return [
                [x, y],
                [x + width - width / 100, y],
                [x + width, y - height / 2],
                [x + width - width / 100, y - height],
                [x, y - height],
            ];
        } else {
            return [
                [x, y - height / 2],
                [x + width / 100, y],
                [x + width, y],
                [x + width, y - height],
                [x + width / 100, y - height],
            ];
        }
    } else if (feature === 'amplicon_feature') {
        return [
            [x, y],
            [x + width, y],
            [x + width, y - height],
            [x, y - height],
        ];
    }
    else {
        return null;
    }
}

/**
 * Renders the gene/amplicon features shape based on below information of geneFeatureData and ampliconData
 * @param params - Echarts arg
 * @param api - Echarts arg
 */
function renderGeneFeatures(params, api) {
    var points, shape, rotateAngle;
    var start, end, height, width, x, y;
    var categoryIndex = params.dataIndex;
    var feature = geneFeatureData[categoryIndex];
    start = api.coord([feature.value.start, categoryIndex]);
    if (categoryIndex === 0) {
        yStart = start[1];
    }
    end = api.coord([feature.value.end, categoryIndex]);
    height = geneFeatureProperties["rec_items_height"];
    width = end[0] - start[0];
    x = start[0];
    y = yStart - height / 2 - feature.value.level;
    points = renderPoints(x, y, width, height, feature.value.strand, feature.value.type);
    if (feature.value.type === 'gene_feature') {
        rotateAngle = (feature.value.strand === 1) ? 0.7 : -0.7;
        shape = echarts.graphic.clipPointsByRect(points, {
            x: params.coordSys.x,
            y: params.coordSys.y,
            width: params.coordSys.width,
            height: params.coordSys.height,
        });
        return {
            type: "polygon",
            shape: {
                points: shape,
            },
            style: api.style({}),
            textContent: {
                type: "text",
                invisible: !showGeneLabel,
                style: {
                    text: geneFeatureData[categoryIndex].name,
                    fill: geneFeatureData[categoryIndex].itemStyle.color,
                    fontStyle: "normal",
                    fontSize: 10,
                    fontWeight: "bolder",
                },
            },
            textConfig: {
                position: "top",
                distance: 20,
                rotation: rotateAngle,
                local: true,
            },
        };
    } else if (feature.value.type === 'amplicon_feature') {
        shape = echarts.graphic.clipPointsByRect(points, {
            x: params.coordSys.x,
            y: params.coordSys.y,
            width: params.coordSys.width,
            height: params.coordSys.height,
        });
        return {
            type: "polygon",
            shape: {
                points: shape,
            },
            style: api.style(),
            textContent: {},
            textConfig: {},
        };
    }
    else {
        return null;
    }
}

/**
 * Define options for gene features charts
 * @param {number} index - gene feature is displayed in the last index of grid
 * @returns {Array<Dict[]>}
 */
function getGeneFeatureSeries(index) {
    var featureSeries = [];
    if (amplicon === 'True' || geneFeature === 'True'){
        featureSeries.push({
            type: "custom",
            xAxisIndex: index,
            yAxisIndex: index,
            renderItem: renderGeneFeatures,
            labelLayout: {
                hideOverlap: false,
            },
            data: geneFeatureData,
            tooltip: {
                trigger: "item",
                enterable: true,
                appendToBody: true,
                renderMode: "html",
                borderRadius: 6,
                borderWidth: 2,
                showContent: "true",
                position: 'top',
                textStyle: {
                    fontSize: 15,
                    fontWeight: "bolder",
                },
                formatter: function (params) {
                    var output = "";
                    rows = [
                        [
                            "Range",
                            params.value.start.toLocaleString() + " - " + params.value.end.toLocaleString(),
                        ],
                        ["Length", (params.value.end - params.value.start + 1).toLocaleString()],
                        ["Strand", params.value.strand],
                    ];
                    output += toTableHtml([params.name, ""], rows, "table small");
                    return output;
                },
            },
        });
    }
    return featureSeries;
}

/**
 * Define options for x Axis
 * @param {Array<string>} samples - An array of samples name
 * @param {number} xAxisMax - Max value is set for x Axis
 * @returns {Array<Dict[]>}
 */
function getXAxes(samples, xAxisMax) {
    var axes = [];
    for (var [i, sample] of samples.entries()) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                interval: "auto",
            },
        });
    }
    if (amplicon === "True" || geneFeature === "True"){
        axes.push({
            type: "value",
            gridIndex: samples.length,
            min: 1,
            max: xAxisMax,
            axisLabel: {
                interval: "auto",
            },
        });
    }
    return axes;
}

/**
 * Define options for Y axis
 * @param {Array<string>} samples - An array of samples name
 * @param {string} scaleType - scale for Y Axis, either value or log
 * @param {number} yMax - max value is set for y Axis
 * @returns {Array<Dict[]>}
 */
function getYAxes(samples, scaleType, yMax) {
    var axes = [];
    for (var [i, sample] of samples.entries()) {
        axes.push({
            type: scaleType,
            gridIndex: i,
            name: sample,
            nameTextStyle: {
                fontStyle: "normal",
                fontWeight: "bolder",
            },
            nameLocation: "end",
            min: 1,
            max: yMax,
            minorSplitLine: {
                show: true,
            },
        });
    }
    if (amplicon === "True" || geneFeature === "True"){
        axes.push({
            max: geneFeatureProperties["max_grid_height"],
            gridIndex: samples.length,
            show: false,
        });
    }
    return axes;
}

/**
 * Get dataset information
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<number>} positions - an array of genome positions which represent in X Axis
 * @returns {Array<Dict[]>}
 */
function getDatasets(depths, positionArray) {
    var datasets = [];
    for (var [i, depthArray] of depths.entries()) {
        datasets.push({
            dimensions: [
                {name: "depth", type: "float"},
                {name: "position", type: "int"},
            ],
            source: {
                position: positionArray,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples name
 * @returns {Array<Dict[]>}
 */
function getDepthSeries(samples) {
    var depthSeries = [];
    for (var [i, sample] of samples.entries()) {
        depthSeries.push({
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
        });
    }
    return depthSeries;
}

/**
 * Define options for amplicon depth coverage bars
 * @param {Array<string>} samples - An array of samples name
 * @returns {Array<Dict[]>}
 */
function getAmpliconDepthSeries(samples) {
    var ampliconDepthSeries = [];
    for (var [i, sample] of samples.entries()) {
        ampliconDepthSeries.push({
            type: "custom",
            xAxisIndex: i,
            yAxisIndex: i,
            renderItem: function (params, api) {
                var start = api.coord([api.value(0), api.value(2)]);
                var end = api.coord([api.value(1), 1]);
                var rectShape = echarts.graphic.clipRectByRect(
                    {
                        x: start[0],
                        y: start[1],
                        width: end[0] - start[0],
                        height: end[1] - start[1]
                    },
                    {
                        x: params.coordSys.x,
                        y: params.coordSys.y,
                        width: params.coordSys.width,
                        height: params.coordSys.height
                    }
                );
                return rectShape && {
                    type: "rect",
                    shape: rectShape,
                    style: api.style({})
                }
            },
            label: {
                show: false,
                position: "top",
                distance: 25,
                rotate:60
            },
            labelLayout: {
                hideOverlap: false
            },
            encode: {
                x: [0, 1],
                y: 2,
            },
            data: ampliconData[sample],
        })
    }
    return ampliconDepthSeries;
}

/**
 * Define options for variant bar charts
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @param {Array<Array<number>>} depths - Array of depths
 * @returns {Array<Dict[]>}
 */
function getVariantsSeries(variants, depths) {
    var variantsSeries = [];
    for (var [i, varMap] of variants.entries()) {
        (function (i, varMap) {
            variantsSeries.push({
                type: "bar",
                xAxisIndex: i,
                yAxisIndex: i,
                data: Object.values(varMap).map((x) => [parseInt(x['POS']), depths[i][x['POS']]]),
                barWidth: 2,
                itemStyle: {
                    color: function (params) {
                        var pos = params.data[0];
                        var nt = window.refSeq[pos-1]
                        if (ntColor.hasOwnProperty(nt)) {
                            return ntColor[nt];
                        }
                        return "#333";
                    },
                },

            });
        })(i, varMap);
    };
    return variantsSeries;
}

/**
 *  Define grid for whole charts
 * @param {Array<string>} samples - An array of samples
 * @returns {Array<Dict[]>}
 */
function getGrids(samples) {
    var n = samples.length + 1;
    var lastHeight;
    var lastTop;
    var grids = Object.keys(samples).map(function (sample) {
        lastHeight = (1 / n) * 100 - 6;
        if (n == 2) {
            // Only 1 sample (1 sample + gene feature plot)
            lastHeight = 70;
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
        var padTop = 4;
        lastTop = (idx / n) * 100 + padTop;
        grid.top = (idx / n) * 100 + padTop + "%";
        grid.left = "8%";
        grid.right = "8%";
    });
    if (amplicon === "True" || geneFeature === "True"){
        grids.push({
            show: true,
            height: geneFeatureProperties["grid_height"],
            top: lastHeight + lastTop + 3 + "%",
            left: "8%",
            right: "8%",
        });
    }
    return grids;
}

/**
 * Define options for tooltips
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @returns {Array<Dict[]>}
 */
function getTooltips(samples, depths, variants) {
    return [
        {
            trigger: "axis",
            enterable: true,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            position: function (pos, params, dom, rect, size) {
                // tooltip will be fixed on the right if mouse hovering on the left,
                // and on the left if hovering on the right.
                var obj = {top: 5};
                obj[['left', 'right'][+(pos[0] < size.viewSize[0] / 2)]] = 5;
                return obj;
            },
            formatter: function (params) {
                var output = "";
                var param = params[0];
                var i = param.axisIndex;
                if (i > samples.length) {
                    return output;
                }
                var sample = samples[i];
                var position = param.data[1];
                var depth = param.data[0];
                var zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
                var zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
                var meanCov = meanCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                var medianCov = medianCoverage(depths, zoomStart, zoomEnd, i).toFixed(2);
                var genomeCov = genomeCoverage(depths, zoomStart, zoomEnd, i, 10).toFixed(2);
                output += "<h5>" + sample + "</h5>";
                var rows = [
                    ["Position", position.toLocaleString()],
                    ["Depth", depth.toLocaleString()],
                ];
                if (params.length > 1) {
                    Object.values(variants[i]).forEach(values => {
                        if (values['POS'] === position) {
                            for (const [key, value] of Object.entries(values)) {
                                if (key !== 'POS' && key !== 'sample') {
                                    rows.push(
                                        ...[[key, value]]
                                    )
                                }

                            }
                        }

                    })
                } else {
                    rows.push(["Sequence", window.refSeq[position - 1]]);
                }
                output += toTableHtml(["Position Info", ""], rows, "table small");
                rows = [
                    [
                        "Range",
                        zoomStart.toLocaleString() + " - " + zoomEnd.toLocaleString(),
                    ],
                    ["Mean Coverage", meanCov + "X"],
                    ["Median Coverage", medianCov + "X"],
                    ["Genome Coverage ( >= 10x)", genomeCov + "%"],
                ];
                output += toTableHtml(["Coverage View Stats", ""], rows, "table small");

                return output;
            },
        },
    ];
}

/**
 * Define all options for coverage chart
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Dict[string, Dict[]]} variants - The dict of variants data
 * @returns {Dict[]}
 */
function getCoverageChartOption(samples, depths, variants) {
    gridLength = samples.length;
    var chartOptions = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, refSeqLength),
        yAxis: getYAxes(samples, "log", maxDepth),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths),
            ...getAmpliconDepthSeries(samples),
            ...getGeneFeatureSeries(gridLength)
        ],
        tooltip: getTooltips(samples, depths, variants),
        toolbox: {
            show: "true",
            feature: {
                dataView: {
                    readOnly: false,
                },
                restore: {},
                saveAsImage: {
                    name: "wgscovplot",
                },
            },
        },
        dataZoom: [
            {
                type: "inside",
                filterMode: "none",
                xAxisIndex: [...Array(gridLength + 1).keys()],
                zoomLock: false,
            },
            {
                show: true,
                filterMode: "none",
                xAxisIndex: [...Array(gridLength + 1).keys()],
                type: "slider",
                zoomLock: false,
            },
        ],
        grid: getGrids(samples)
    };
    return chartOptions;
}

/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are still reserved
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
    var geneFeatureHeight;
    samples.forEach(sample => {
        depths.push(window.depths[sample]);
        variants.push(window.variants[sample]);
    })
    /* Need to handle the case when user removes all samples on the chart and add new samples to be displayed
    Variables lenCheck1, lenCheck2 and closures lastYAxisScale, lastYAxisMax are used to keep track the scale settings of a last sample on the chart before it was removed
    lenCheck1, lenCheck1 depends on whether the last subplot is used for plotting gene/amplicon features or not
     */
    var chartOption = chart.getOption();
    var lenCheck1 = (amplicon === 'True' || geneFeature === 'True') ? 2 : 1
    var lenCheck2 = (amplicon === 'True' || geneFeature === 'True') ? 1 : 0
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
    //Left/right margin
    if (chartOption.grid.length === 0){ // There is no subplot and gene/amplicon feature
        leftMargin = 8
        rightMargin = 8
        document.getElementById("chart-left-input").value = _.toInteger(leftMargin);
        document.getElementById("chart-left-output").value = _.toInteger(leftMargin) + "%";
        document.getElementById("chart-right-input").value = _.toInteger(rightMargin);
        document.getElementById("chart-right-output").value = _.toInteger(rightMargin) + "%";
    }
    else{
        leftMargin = _.toInteger(chartOption.grid[0].left.replace("%", ""));
        rightMargin = _.toInteger(chartOption.grid[0].right.replace("%", ""));
    }
    // Get Current Data Zoom
    var zoomStart = Math.floor(chart.getOption().dataZoom[0].startValue);
    var zoomEnd = Math.floor(chart.getOption().dataZoom[0].endValue);
    // The current chart is not disposed so notMerge must be set true
    chart.setOption(option = getCoverageChartOption(samples, depths, variants), notMerge = true);
    setScale(scaleType, max);
    setDataZoom(zoomStart, zoomEnd);
    updateChartLeftMargin(leftMargin);
    updateChartRightMargin(rightMargin);
    updateControlMenu();
}

/**
 * When the chart is firstly initialized, the samples to plot on control menu is updated with first 3 samples
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
 * Get the list of current samples in the Samples to plot menu
 * @returns {Array<string>>} An array of samples name
 */
function getCurrentSamples(chartOption) {
    var samples = [];
    // If the chart is not initalizaed yet, get 3 first samples from window.samples
    if (chartOption === undefined || chartOption === null){
        samples = _.slice(window.samples, 0, defaultNumChart);
    } else{
        var selectData = $("#selectedsamples").select2("data");
        for (var [key, entries] of selectData.entries()) {
            samples.push(selectData[key].text);
        }
    }
    return samples;
}

/**
 * Initialize the events handler for chart after it is initialized
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
            _.forEach(seriesOption, function (element) {
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
         * An action updates chart options when number of selected samples changes
         * The chart options such as y Axis scale, yMax, DataZoom are still reserved (users's settings are respected)
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
            seriesOption[seriesOption.length - 1]["renderItem"] = renderGeneFeatures; // Re-update Gene Feature Chart Only
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
            chart.setOption({
                dataZoom: [
                    {
                        type: "inside",
                        filterMode: "none",
                        xAxisIndex: [...Array(gridLength + 1).keys()],
                    },
                    {
                        type: "slider",
                        show: isChecked,
                        filterMode: "none",
                        xAxisIndex: isChecked ? [...Array(gridLength + 1).keys()] : null,
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
        chart.setOption(option = getCoverageChartOption(plotSamples, plotDepths, plotVariants));
    } else {
        var renderEnv = document.getElementById("render-env").value;
        var isChecked = document.getElementById("toggle-darkmode").checked;
        var mode = isChecked ? "dark" : "white";
        var gridOption = chart.getOption().grid;
        var scaleType = chart.getOption().yAxis[0].type;
        var max = chart.getOption().yAxis[0].max;
        var dataZoomOption = chart.getOption().dataZoom;
        echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        chart = echarts.init($chart, mode, {renderer: renderEnv});
        chart.setOption(option = getCoverageChartOption(plotSamples, plotDepths, plotVariants));
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
                _.toInteger(gridOption[i - 1]["top"].replace("%", "")) +
                _.toInteger(gridOption[i - 1]["height"].replace("%", "")) +
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
    _.forEach(gridOption, function (element) {
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
    _.forEach(gridOption, function (element) {
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
                _.toInteger(gridOption[i - 1]["top"].replace("%", "")) +
                _.toInteger(gridOption[i - 1]["height"].replace("%", "")) +
                _.toInteger(val) +
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
        _.forEach(yAxisOption, function (element) {
            if (element.gridIndex < len) {
                element.type = scale;
                element.min = 0;
                element.max = yMax;
            }
        });
    } else {
        _.forEach(yAxisOption, function (element) {
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
    _.forEach(yAxisOption, function (element) {
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
        var height1 = _.replace(gridOption[0].height, "%", "");
        var top = _.replace(gridOption[0].top, "%", "");
        document.getElementById("chart-height-input").value = _.toInteger(height1);
        document.getElementById("chart-height-output").value = _.toInteger(height1) + "%";
        document.getElementById("chart-top-input").value = _.toInteger(top);
        document.getElementById("chart-top-output").value = _.toInteger(top) + "%";
        if (amplicon === 'True' || geneFeature === 'True') {
            var height2 = _.replace(gridOption[gridOption.length - 1].height, "%", "");
            document.getElementById("genefeature-height-input").value = _.toInteger(height2);
            document.getElementById("genefeature-height-output").value = _.toInteger(height2) + "%";
        }
    }
}

/**
 * Write tooltip information to HTML table
 * @param {string} headers - Header of table
 * @param {Array<string>} rows - Rows of table
 * @param {string} classes - Classes defined for table
 * @returns {string}
 */
function toTableHtml(headers, rows, classes) {
    var classTable = _.defaultTo(classes, "table");
    var out = '<table class="' + classTable + '"><thead>';
    out += _.join(
        _.map(headers, function (x) {
            return "<strong>" + x + "</strong>";
        }),
        ""
    );
    out += "</thead><tbody>";
    out += _.join(
        _.map(rows, function (xs) {
            return (
                "<tr>" +
                _.join(
                    _.map(xs, function (x, i) {
                        return "<td " + (i === 0 ? 'scope="row"' : "") + ">" + x + "</td>";
                    }),
                    ""
                ) +
                "</tr>"
            );
        }),
        ""
    );
    out += "</tbody></table>";
    return out;
}

/**
 * Calculate mean coverage for a genome region
 * @param {Array<number>} depths - depths array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @returns {number}
 */
function meanCoverage(depths, start, end, gridIndex) {
    var subArray = _.slice(depths[gridIndex], start - 1, end);
    return _.mean(subArray);
}

/**
 *  Calculate genome coverage depth acorrding to threshod low
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @param {number} low - the threshold that wants to set
 * @returns {number}
 */
function genomeCoverage(depths, start, end, gridIndex, low) {
    var subArray = _.slice(depths[gridIndex], start - 1, end);
    var filetedArray = _.filter(subArray, function (x) {
        return x >= low;
    });
    return (filetedArray.length / (end - start + 1)) * 100;
}

/**
 * Calculate median of an array
 * @param {Array<number>} arr
 * @returns {number}
 */
function median(arr) {
    var sortedArr = _.sortBy(arr)
    var half = Math.floor(sortedArr.length / 2);
    if (sortedArr.length % 2) return sortedArr[half];
    else return (sortedArr[half - 1] + sortedArr[half]) / 2.0;
}

/**
 * Calculate median coverage for a genome region
 * @param {Array<number>} depths - depth array
 * @param {number} start - start position
 * @param {number} end - end position
 * @param {number} gridIndex - grid index of sample
 * @returns {number}
 */
function medianCoverage(depths, start, end, gridIndex) {
    var subArray = _.slice(depths[gridIndex], start - 1, end);
    return median(subArray);
}
