/**
 * The function returns the points of gene features shape
 * @param x - x-axis coordinate
 * @param y  - x-axis coordinate
 * @param width - width of shape
 * @param height - height of shape
 * @param strand - strand of feature
 * @param feature - gene feature or amplicon feature
 * @returns { Coordinate of 5 points for gene feature or 4 points for amplicon feature }
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
            ]
        } else {
            return [
                [x, y - height / 2],
                [x + width / 100, y],
                [x + width, y],
                [x + width, y - height],
                [x + width / 100, y - height],
            ]
        }
    } else if (feature === 'amplicon_feature') {
        return [
            [x, y],
            [x + width, y],
            [x + width, y - height],
            [x, y - height],
        ]
    }
    else {
        return null
    }
}

/**
 * The function renders the gene features shape based on below information of gene_feature
 * @param params - Echarts arg
 * @param api - Echarts arg
 * @returns { shape of gene feature or amplicon feature based on points returned from renderPoints function }
 */
function renderGeneFeatures(params, api) {
    var points, shape, rotate_angle;
    var start, end, height, width, x, y;
    //var categoryIndex = api.value(0);
    var categoryIndex = params.dataIndex;
    var feature = gene_feature[categoryIndex]
    start = api.coord([feature.value.start, categoryIndex]);
    if (categoryIndex === 0) {
        y_start = start[1];
    }
    end = api.coord([feature.value.end, categoryIndex]);
    height = gene_feature_properties["rec_items_height"];
    width = end[0] - start[0];
    x = start[0];
    y = y_start - height / 2 - feature.value.level;
    if (feature.value.type === 'gene_feature') {
        points = renderPoints(x, y, width, height, feature.value.strand, 'gene_feature');
        if (feature.value.strand === 1) { // Plus Strand
            rotate_angle = 0.7;
        } else { // Minus Strand
            rotate_angle = -0.7;
        }
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
                invisible: invisible_gene_label,
                style: {
                    text: gene_feature[categoryIndex].name,
                    fill: gene_feature[categoryIndex].itemStyle.color,
                    fontStyle: "normal",
                    fontSize: 10,
                    fontWeight: "bolder",
                },
            },
            textConfig: {
                position: "top",
                distance: 20,
                rotation: rotate_angle,
                local: true,
            },
        };
    } else if (feature.value.type === 'amplicon_feature') {
        points = renderPoints(x, y, width, height, feature.value.strand, 'amplicon_feature');
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
        return null
    }
}

/**
 * The function returns options for gene features charts
 * @param index - gene feature is displayed in the last index of grid
 * @returns {*[]}
 */
function getGeneFeatureSeries(index) {
    var feature_series = [];
    feature_series.push({
        type: "custom",
        xAxisIndex: index,
        yAxisIndex: index,
        renderItem: renderGeneFeatures,
        labelLayout: {
            hideOverlap: false,
        },
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
    return feature_series;
}

/**
 * The functions returns options for X axis
 * @param samples - the list samples name
 * @param ref_len -the length of reference sequence
 * @returns {*[]}
 */
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

/**
 * The functions returns options for Y axis
 * @param samples - the list of samples name
 * @param scaletype - scale for Y Axis: either linear or log scale
 * @param ymax - max value of y Axis
 * @returns {*[]}
 */
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
        max: gene_feature_properties["max_grid_height"],
        gridIndex: samples.length,
        show: false,
    });
    return axes;
}

/**
 * The function returns dataset information
 * @param depths - the dict of depth data
 * @param positions - an array of genome positions which represent in X Axis
 * @returns {*[]}
 */
function getDatasets(depths, positions) {
    var datasets = [];
    for (var [i, depthArray] of depths.entries()) {
        datasets.push({
            dimensions: [
                {name: "depth", type: "float"},
                {name: "position", type: "int"},
            ],
            source: {
                position: positions,
                depth: depthArray,
            },
        });
    }
    return datasets;
}

/**
 * The function returns options for depth coverage charts
 * @param samples - the list of samples name
 * @returns []
 */
function getDepthSeries(samples) {
    var series = [];
    for (var [i, sample] of samples.entries()) {
        series.push({
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
    return series;
}

/**
 * The function returns options for  amplicon depth coverage bars
 * @param samples - the list of samples name
 * @returns []
 */
function getAmpliconDepthSeries(samples) {
    var series = []
    for (var [i, sample] of samples.entries()) {
        series.push({
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
            data: amplicon_data[sample],
        })
    }
    return series
}

/**
 * The function returns options for variant charts
 * @param variants - the dict of variants data
 * @param depths - the dict of depths data
 * @returns {*[]}
 */
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
                        }
                        return "#333";
                    },
                },
            });
        })(i, varMap);
    }
    return series;
}

/**
 *  The function defines grid for whole charts
 * @param samples - the list of samples name
 * @returns []
 */
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
        var padTop = 4;
        last_top = (idx / n) * 100 + padTop;
        grid.top = (idx / n) * 100 + padTop + "%";
        grid.left = "5%";
        grid.right = "5%";
    });
    grids.push({
        show: true,
        height: gene_feature_properties["grid_height"],
        top: last_height + last_top + 3 + "%",
        left: "5%",
        right: "5%",
    });
    return grids;
}

/**
 * The function defines options for tooltips
 * @param samples - the list of sample names
 * @param depths - the dict of depths data
 * @param variants - the dict of variants data
 * @returns {[{renderMode: string, formatter: ((function(*): (string))|*), enterable: boolean, appendToBody: boolean, showContent: boolean, trigger: string}]}
 */
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
                if (i > samples.length) {
                    return output;
                }
                var sample = samples[i];
                var position = param.data[1];
                var depth = param.data[0];
                var start_pos = Math.floor(chart.getOption().dataZoom[0].startValue);
                var end_pos = Math.floor(chart.getOption().dataZoom[0].endValue);
                var mean_cov = meanCoverage(depths, start_pos, end_pos, i).toFixed(2);
                var median_cov = medianCoverage(depths, start_pos, end_pos, i).toFixed(2);
                var genome_cov = genomeCoverage(depths, start_pos, end_pos, i, 10).toFixed(2);
                output += "<h5>" + sample + "</h5>";
                var rows = [
                    ["Position", position.toLocaleString()],
                    ["Depth", depth.toLocaleString()],
                ];
                if (variants[i].hasOwnProperty(position)) {
                    rows.push(
                        ...[
                            [
                                "Ref",
                                ref_seq.substring(
                                    position - 1,
                                    position - 1 + variants[i][position].length
                                ),
                            ],
                            ["Variant", variants[i][position]],
                        ]
                    );
                } else {
                    rows.push(["Sequence", ref_seq[position - 1]]);
                }
                output += toTableHtml(["Position Info", ""], rows, "table small");
                rows = [
                    [
                        "Range",
                        start_pos.toLocaleString() + " - " + end_pos.toLocaleString(),
                    ],
                    ["Mean Coverage", mean_cov + "X"],
                    ["Median Coverage", median_cov + "X"],
                    ["Genome Coverage ( >= 10x)", genome_cov + "%"],
                ];
                output += toTableHtml(["Coverage View Stats", ""], rows, "table small");
                return output;
            },
        },
    ];
}

/**
 * The function builds all options for coverage chart
 * @param samples - the list of sample names
 * @param depths - the dict of depths data
 * @param variants - the dict of variants data
 * @returns {yAxis: *[], xAxis: *[], series: *[], grid: [], tooltip: {renderMode: string, formatter: ((function(*): string)|*), enterable: boolean, appendToBody: boolean, showContent: boolean, trigger: string}[], toolbox: {feature: {saveAsImage: {name: string}, restore: {}, dataView: {readOnly: boolean}}, show: string}, dataZoom: [{filterMode: string, xAxisIndex: number[], type: string, zoomLock: boolean}, {filterMode: string, xAxisIndex: number[], show: boolean, type: string, zoomLock: boolean}], title: {}, dataset: *[]}
 */
function getCoverageChartOption(samples, depths, variants) {
    grid_length = samples.length;
    var options = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, ref_len),
        yAxis: getYAxes(samples, "log", max_depth),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths),
            ...getAmpliconDepthSeries(samples),
            ...getGeneFeatureSeries(grid_length)
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
                    name: "Coverage_Plot",
                },
            },
        },
        dataZoom: [
            {
                type: "inside",
                filterMode: "none",
                xAxisIndex: [...Array(grid_length + 1).keys()],
                zoomLock: false,
            },
            {
                show: true,
                filterMode: "none",
                xAxisIndex: [...Array(grid_length + 1).keys()],
                type: "slider",
                zoomLock: false,
            },
        ],
        grid: getGrids(samples),

    };
    return options;
}

/**
 * The function updates options for coverage charts.
 * Whenever the selected samples changed, the whole chart will be re-rendered but not dispose the chart
 * @param samples -the list samples name
 */
function updateCoverageChartOption(samples) {
    var depths = [];
    var variants = [];
    samples.forEach(sample => {
        depths.push(window.depths[sample]);
        variants.push(window.variants[sample]);
    })
    // The current chart is not disposed so norMerge must be set true
    chart.setOption(option = getCoverageChartOption(samples, depths, variants), notMerge = true);
    updateControlMenu();
}

/**
 * The function returns the list of first 'default_num_chart = 3' samples when the chart is initialized
 * @param samples - the list samples name
 */
function selectDefaultSamples(updated_samples) {
    // Set default samples display
    $("#selectedsamples").select2();
    $("#selectedsamples").val(updated_samples);
    $("#selectedsamples").trigger('change')
}

/**
 * The function to get the current samples of multiple select of select2
 * @returns {*[]}
 */
function selectCurrentSamples() {
    var selectData = $("#selectedsamples").select2("data");
    var current_samples = [];
    for (var [key, entries] of selectData.entries()) {
        current_samples.push(selectData[key].text);
    }
    return current_samples
}

/**
 * A document is ready for fire action with Jquery
 */
$(document).ready(function () {

    $("#selectedsamples").select2({
        tags: true,
    });

    /**
     * Jquery function to make the list of samples is not forced in alphabetical order
     */
    $("#selectedsamples").on("select2:select", function (evt) {
        var element = evt.params.data.element;
        var $element = $(element);
        $element.detach();
        $(this).append($element);
        $(this).trigger("change");
    });

    /**
     * Jquery function update chart options when number of selected samples changes
     */
    $("#selectedsamples").on("change", function () {
        updateCoverageChartOption(selectCurrentSamples());
    });

    /**
     * Toggle dark mode, the entire chart will be disposed and re-initialized
     * And the current selected samples of select2 will re-rendered
     */
    $("#toggle-darkmode").change(function () {
        var render_env = "canvas";
        if (document.getElementById("render-env"))
            render_env = document.getElementById("render-env").value;
        var current_samples = selectCurrentSamples();
        var current_depths = []
        var current_variants = []
        current_samples.forEach(sample => {
            current_depths.push(window.depths[sample]);
            current_variants.push(window.variants[sample]);
        })

        echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        if ($(this).prop("checked")) {
            chart = echarts.init($chart, "dark", {renderer: render_env});
        } else {
            chart = echarts.init($chart, "white", {renderer: render_env});
        }
        chart.setOption(option = getCoverageChartOption(current_samples, current_depths, current_variants));
        chartDispatchDataZoomAction();
        updateControlMenu();
    });

    /**
     * Toggle to Amplicon Depth Label
     */
    $("#toggle-amplicon-depthlabel").change(function () {
        var series_option = chart.getOption().series;
        var isChecked = $(this).prop("checked")
        _.forEach(series_option, function (element) {
            if (element.type === 'custom') {
                element.label.show = isChecked
            }
        })
        chart.setOption({series: [...series_option]});
    });

    /**
     * Toggle to show gene label or not
     */
    $("#toggle-genelabel").change(function () {
        var series_option = chart.getOption().series;
        invisible_gene_label = !($(this).prop("checked"))
        series_option[series_option.length - 1]["renderItem"] = renderGeneFeatures; // Re-update Gene Feature Chart Only
        chart.setOption({series: [...series_option]});
    });

    /**
     * Toggle tooltip for coverage chart, this action does not affect gene feature chart
     */
    $("#toggle-tooltip").change(function () {
        var isChecked = $(this).prop("checked")
        chart.setOption({
            tooltip: {showContent: isChecked},
        });
    });

    /**
     * Toggle slider zoom
     */
    $("#toggle-slider").change(function () {
        var isChecked = $(this).prop("checked")
        chart.setOption({
            dataZoom: [
                {
                    type: "inside",
                    filterMode: "none",
                    xAxisIndex: [...Array(grid_length + 1).keys()],
                },
                {
                    show: isChecked,
                    filterMode: "none",
                    xAxisIndex: isChecked ? [...Array(grid_length + 1).keys()] : null,
                    type: "slider",
                },
            ],
        })
    });
});

/**
 * Select svg or canvas as rendered environment, the entire chart will be disposed and re-initialized
 * And the current selected samples of select2 will re-rendered
 */
function selectRenderEnv() {
    var render_env = document.getElementById("render-env").value;
    var isChecked = document.getElementById("toggle-darkmode").checked;
    var current_samples = selectCurrentSamples();
    var current_depths = []
    var current_variants = []
    current_samples.forEach(sample => {
        current_depths.push(window.depths[sample]);
        current_variants.push(window.variants[sample]);
    })
    var mode = "white";
    if (isChecked) mode = "dark";
    else mode = "white";
    echarts.dispose(chart); // destroy chart instance and re-init chart
    $chart = document.getElementById("chart");
    if (render_env === "canvas") {
        chart = echarts.init($chart, mode, {renderer: "canvas"});
    } else {
        chart = echarts.init($chart, mode, {renderer: "svg"});
    }
    chart.setOption(option = getCoverageChartOption(current_samples, current_depths, current_variants));
    chartDispatchDataZoomAction();
    updateControlMenu();
}

/**
 * Adjust chart height
 * @param val - Subplots height percent value
 */
function updateChartHeight(val) {
    document.getElementById("chart-height-output").value = val + "%";
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
    chart.setOption({grid: grid_option});
}

/**
 * Adjust left margin of chart
 * @param val - Subplots left margin percent value
 */
function updateChartLeft(val) {
    document.getElementById("chart-left-output").value = val + "%";
    var grid_option = chart.getOption().grid;
    _.forEach(grid_option, function (element) {
        element.left = val + "%";
    });
    chart.setOption({grid: grid_option});
}

/**
 * Adjust right margin of chart
 * @param val - Subplots right margin percent value
 */
function updateChartRight(val) {
    document.getElementById("chart-right-output").value = val + "%";
    var grid_option = chart.getOption().grid;
    _.forEach(grid_option, function (element) {
        element.right = val + "%";
    });
    chart.setOption({grid: grid_option});
}

/**
 * Adjust top margin of chart
 * @param val - Subplots top margin percent value
 */
function updateChartTop(val) {
    document.getElementById("chart-top-output").value = val + "%";
    var grid_option = chart.getOption().grid;
    for (var i = 0; i < grid_option.length; i++) {
        if (i === 0) {
            grid_option[i]["top"] = val + "%";
        } else {
            grid_option[i]["top"] =
                _.toInteger(grid_option[i - 1]["top"].replace("%", "")) +
                _.toInteger(grid_option[i - 1]["height"].replace("%", "")) +
                _.toInteger(val) +
                "%";
        }
    }
    chart.setOption({grid: grid_option});
}



/**
 * Adjust the height of gene feature charts
 * @param val - Subplots gene feature height percent value
 */
function updateGeneFeatureHeight(val) {
    document.getElementById("genefeature-height-output").value = val + "%";
    var grid_option = chart.getOption().grid;
    gene_features_grid_index = grid_option.length - 1;
    grid_option[gene_features_grid_index]["height"] = val + "%";
    chart.setOption({grid: grid_option});
}

/**
 * Select log or linear scale for the chart
 */
function setScale() {
    var scale_type = document.getElementById("scale").value;
    var yaxis_option = chart.getOption().yAxis;
    if (scale_type === "value") {
        _.forEach(yaxis_option, function (element) {
            if (element.gridIndex < yaxis_option.length - 1) {
                element.type = scale_type;
                element.min = 0;
                element.max = max_depth;
            }
        });
        chart.setOption({yAxis: yaxis_option});
        document.getElementById("ymax").value = max_depth; // update control menu
    } else {
        _.forEach(yaxis_option, function (element) {
            if (element.gridIndex < yaxis_option.length - 1) {
                element.type = scale_type;
                element.min = 1;
                element.max = max_depth;
            }
        });
        chart.setOption({yAxis: yaxis_option});
        document.getElementById("ymax").value = max_depth; // update control menu
    }
}

/**
 * Set YMax for Y Axis
 */
function setYMax() {
    var ymax = document.getElementById("ymax").value;
    var yaxis_option = chart.getOption().yAxis;
    _.forEach(yaxis_option, function (element) {
        if (element.gridIndex < yaxis_option.length - 1) {
            element.max = ymax;
        }
    });
    chart.setOption({yAxis: yaxis_option});
}

/**
 * Set range of gene feature to zoom in
 */
function setDataZoom() {
    var start = document.getElementById("start-pos").value;
    var end = document.getElementById("end-pos").value;
    chart.dispatchAction({
        type: "dataZoom",
        startValue: start,
        endValue: end,
    });
}

/**
 * Reset range of gene feature to default [1, ref_len]
 */
function resetDataZoom() {
    document.getElementById("start-pos").value = 1;
    document.getElementById("end-pos").value = ref_len;
    chart.dispatchAction({
        type: "dataZoom",
        start: 0,
        end: 100,
    });
}


/**
 * The Control Menu is updated when render env, toggle dark mode or the selected samples changed
 */
function updateControlMenu() {
    var grid_option = chart.getOption().grid;
    var height1 = _.replace(grid_option[0].height, "%", "");
    var height2 = _.replace(grid_option[grid_option.length - 1].height, "%", "");
    var left = _.replace(grid_option[0].left, "%", "");
    var right = _.replace(grid_option[0].right, "%", "");
    var top = _.replace(grid_option[0].top, "%", "");
    // Reset Chart Properties and Gene Features Menu
    document.getElementById("chart-height-input").value = _.toInteger(height1);
    document.getElementById("chart-height-output").value = _.toInteger(height1) + "%";
    document.getElementById("chart-left-input").value = _.toInteger(left);
    document.getElementById("chart-left-output").value = _.toInteger(left) + "%";
    document.getElementById("chart-right-input").value = _.toInteger(right);
    document.getElementById("chart-right-output").value = _.toInteger(right) + "%";
    document.getElementById("chart-top-input").value = _.toInteger(top);
    document.getElementById("chart-top-output").value = _.toInteger(top) + "%";
    document.getElementById("genefeature-height-input").value = _.toInteger(height2);
    document.getElementById("genefeature-height-output").value = _.toInteger(height2) + "%";
    // Set Axis to Log scale
    document.getElementById("scale").value = "log"
    document.getElementById("ymax").value = max_depth
    // Reset Data Zoom
    resetDataZoom()
}

/**
 * Dispacth click, dbclick event to DataZoom
 * The function will be called when the chart is intialized, toggle dark mode and change render env
 */
function chartDispatchDataZoomAction() {
    chart.on("click", function (params) {
        document.getElementById("start-pos").value = params.value.start;
        document.getElementById("end-pos").value = params.value.end;
        chart.dispatchAction({
            type: "dataZoom",
            startValue: params.value.start,
            endValue: params.value.end,
        });
    });

    chart.on("dblclick", function (params) {
        document.getElementById("start-pos").value = 1;
        document.getElementById("end-pos").value = ref_len;
        chart.dispatchAction({
            type: "dataZoom",
            start: 0,
            end: 100,
        });
    });
}

/**
 * Write tooltip information to HTML table
 * @param headers - Header of table
 * @param rows - Rows of table
 * @param classes - Classes defined for table
 * @returns {string}
 */
function toTableHtml(headers, rows, classes) {
    classes = _.defaultTo(classes, "table");
    var out = '<table class="' + classes + '"><thead>';
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
 * Calculate mean coverage for specific range
 * @param depths - depth array
 * @param start - start position
 * @param end - end position
 * @param gridIndex - grid index of sample
 * @returns {number}
 */
function meanCoverage(depths, start, end, gridIndex) {
    var sub_array = _.slice(depths[gridIndex], start - 1, end);
    return _.mean(sub_array);
}

/**
 *  Calculate genome coverage depth which is >= low (10)
 * @param depths - depth array
 * @param start - start position
 * @param end - end position
 * @param gridIndex - grid index of sample
 * @param low
 * @returns {number}
 */
function genomeCoverage(depths, start, end, gridIndex, low) {
    var sub_array = _.slice(depths[gridIndex], start - 1, end);
    var fileted_array = _.filter(sub_array, function (x) {
        return x >= low;
    });
    return (fileted_array.length / (end - start + 1)) * 100;
}

/**
 * Calculate median of an array
 * @param arr
 * @returns {number}
 */
function median(arr) {
    arr.sort(function (a, b) {
        return a - b;
    });
    var half = Math.floor(arr.length / 2);
    if (arr.length % 2) return arr[half];
    else return (arr[half - 1] + arr[half]) / 2.0;
}

/**
 * Calculate median coverage for specific range
 * @param depths - depth array
 * @param start - start pos
 * @param end - end pos
 * @param gridIndex - grid index of sample
 * @returns {number}
 */
function medianCoverage(depths, start, end, gridIndex) {
    var sub_array = _.slice(depths[gridIndex], start - 1, end);
    return median(sub_array);
}
