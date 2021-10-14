/**
 * Reference Sequence
 */
window.ref_seq = "{{ ref_seq }}";

/**
 * Dict of samples name
 */
window.samples = {{ samples_name }}

/**
 * Dict of depth data
 * {samples_name: depth_data}
 */
window.depths = {{ depth_data }};

/**
 * Dict of variant data
 * {samples_name: variant_data}
 */
window.variants = {{ variant_data }};

/**
 * Dict of gene features
 * [{
 *     'name': 'feature_name'
 *     'value: [index, start_pos, end_pos, level, strand, 'gene_features' or 'amplicon_features']
 *     'itemStyle': {
 *          'color': 'color'
 *     }
 * }]
 */
var gene_feature = {{ gene_feature }}

/**
 * Default properties for gene features chart, to have best view for the whole, user can adjust high, right/left/top margin
 *  'max_grid_height': 80,
 *  'rec_items_height': 12,
 *  'minus_strand_level': 55,
 *  'grid_height': "15%"
 *
 */
var gene_feature_properties = {{ gene_feature_properties }}

/**
 * Whether amplicon plot or normal coverage plot
 * @type {boolean}
 */
var amplicon = "{{amplicon}}"

/**
 * Dict of Amplicon Depth Coverage
 * { sample_name: [
 *      {   value : [start, end, depth, amplicon_name],
 *          itemStyle: { color: 'color'
 *      }
 *    ]
 * }
 *
 */
var amplicon_data = {{ amplicon_data }}
const ref_len = window.ref_seq.length;

/**
 * Variable is used keep track the current number samples in the chart and is used for toggling slider
 * @type {number}
 */
var grid_length;

/**
 * Variable is used for initial y_start for rendering gene features
 * @type {number}
 */
var y_start;

const default_num_chart = 3;

/**
 * An array of positions of reference sequence which represent in X-Axis
 * @type {number[]}
 */
var positions = [...Array(ref_len + 1).keys()];
positions.shift();

/**
 * Variable is used to toggle gene label name
 * @type {boolean}
 */
var invisible_gene_label = false;

/**
 * Dict of color for coloring variant position in the chart
 * @type {A: string, C: string, T: string, G: string}
 */
var ntColor = {
    "A": "#F00",
    "C": "#0F0",
    "G": "#00F",
    "T": "#0FF",
}
var $chart = document.getElementById('chart');
window.$chart = $chart;
var chart = echarts.init($chart, 'white', {renderer: 'canvas'});
window.chart = chart;

/**
 * The function returns the points of gene features shape
 * @param x, y, width, height, strand, feature
 * @returns ((*|number)[]|(*)[]|(*|number)[])[]|(*[]|(number|*)[]|(*|number)[]|number[]|(*|number)[])[]
 */
function renderPoints(x, y, width, height, strand, feature = "gene_feature") {
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
}
/**
 * The function renders the gene features shape based on below information of gene_feature
 * @param params, api
 * @returns Shape of gene features {shape: {points}, textConfig: {distance: number, rotation: number, position: string, local: boolean}, style, textContent: {invisible: boolean, style: {fontSize: number, text, fill: ((function(*): (*|string))|*), fontStyle: string, fontWeight: string}, type: string}, type: string}
 */
function renderGeneFeatures(params, api) {
    var points, shape, rotate_angle;
    var start, end, height, width, x, y;
    var categoryIndex = api.value(0);
    start = api.coord([api.value(1), categoryIndex]);
    if (categoryIndex === 0) {
        y_start = start[1];
    }
    end = api.coord([api.value(2), categoryIndex]);
    height = gene_feature_properties["rec_items_height"];
    width = end[0] - start[0];
    x = start[0];
    y = y_start - height / 2 - api.value(3);
    if (api.value(5) === 'gene_feature') {
        points = renderPoints(x, y, width, height, api.value(4), 'gene_feature');
        if (api.value(4) === 1) { // Plus Strand
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
    } else if (api.value(5) === 'amplicon_feature') {
        points = renderPoints(x, y, width, height, api.value(4), 'amplicon_feature');
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
}

/**
 * The function returns options for gene features charts
 * @param index (the last index grid index array)
 * @returns []
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
                    (params.value[2] - params.value[1] + 1).toLocaleString() +
                    "<br/>" +
                    "Strand: " +
                    params.value[4].toLocaleString();
                return output;
            },
        },
    });
    return feature_series;
}

/**
 * The functions returns options for X axis
 * @param samples, ref_len
 * @returns []
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
 * @param samples, scaletype, ymax
 * @returns []
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
 * @param depths, positions
 * @returns []
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
 * @param samples
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
 * @param samples
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
 * @param variants, depths
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
 * @param samples
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
 * @param samples, depths, variants
 * @returns [{renderMode: string, formatter: ((function(*): (string))|*), enterable: boolean, appendToBody: boolean, showContent: boolean, trigger: string}]
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
 * @returns {yAxis: [], xAxis: [], series: *[], grid: [], tooltip: {renderMode: string, formatter: ((function(*): string)|*), enterable: boolean, appendToBody: boolean, showContent: boolean, trigger: string}[], toolbox: {feature: {saveAsImage: {name: string}, restore: {}, dataView: {readOnly: boolean}}, show: string}, dataZoom: [{filterMode: string, xAxisIndex: number[], type: string, zoomLock: boolean}, {filterMode: string, xAxisIndex: number[], show: boolean, type: string, zoomLock: boolean}], title: {}, dataset: []}
 */
function getCoverageChartOption() {
    var samples = [];
    var depths = [];
    var variants = [];

    for (const [key, entries] of Object.entries(window.samples)) {
        if (key < default_num_chart) {
            samples.push(entries);
            depths.push(window.depths[entries]);
            variants.push(window.variants[entries]);
        }
    }
    grid_length = samples.length;
    var options = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, ref_len),
        yAxis: getYAxes(samples, "log", 100000),
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
    selectDefaultSamples(samples);
    return options;
}

/**
 * The function updates options for coverage charts.
 * Whenever the selected samples changed, the whole chart will be re-rendered
 * @param samples
 */
function updateCoverageChartOption(samples) {
    var depths = [];
    var variants = [];
    for (const selected_samples of samples) {
        depths.push(window.depths[selected_samples]);
        variants.push(window.variants[selected_samples]);
    }
    grid_length = samples.length;
    var options = {
        title: {},
        dataset: getDatasets(depths, positions),
        xAxis: getXAxes(samples, ref_len),
        yAxis: getYAxes(samples, "log", 100000),
        // Render 1. Coverage depth; 2. Variants; 3 Amplicon Bar Plot; 4. Gene Feature
        series: [
            ...getDepthSeries(samples),
            ...getVariantsSeries(variants, depths),
            ...getAmpliconDepthSeries(samples),
            ...getGeneFeatureSeries(grid_length),
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
    updateControlMenu();
}

chart.setOption((option = getCoverageChartOption()));

/**
 * The function returns the list of first 'default_num_chart = 3' samples when the first first initialized
 * @param samples
 */
function selectDefaultSamples(samples) {
    // Set default samples display
    $("#selectedsamples").select2();
    $("#selectedsamples").val(samples);
}

/**
 * Jquery function update chart options when number of selected samples changes
 */
$(document).ready(function () {
    $("#selectedsamples").select2();
    $("#selectedsamples").on("change", function (e) {
        var selectData = $("#selectedsamples").select2("data");
        var samples_list = [];
        for (var [key, entries] of selectData.entries()) {
            samples_list.push(selectData[key].text);
        }
        updateCoverageChartOption(samples_list);
    });
});

/**
 * Jquery function to make the list of samples is not forced in alphabetical order
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
 * Select svg or canvas as rendered environment, the entire chart/env will be re-rendered
 */
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
    echarts.dispose(chart); // destroy chart instance and re-init chart
    $chart = document.getElementById("chart");
    if (render_env === "canvas") {
        chart = echarts.init($chart, mode, {renderer: "canvas"});
    } else {
        chart = echarts.init($chart, mode, {renderer: "svg"});
    }
    chart.setOption(option = getCoverageChartOption());
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateControlMenu();
}

/**
 * Toggle dark mode, the entire chart/env will be re-rendered
 */
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
    echarts.dispose(chart); // destroy chart instance and re-init chart
    $chart = document.getElementById("chart");
    if ($(this).prop("checked")) {
        chart = echarts.init($chart, "dark", {renderer: render_env});
    } else {
        chart = echarts.init($chart, "white", {renderer: render_env});
    }
    chart.setOption(option = getCoverageChartOption());
    selectDefaultSamples(samples);
    chartDispatchAction();
    updateControlMenu();
});

/**
 * Adjust chart height
 * @param val
 */
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
    chart.setOption({grid: grid_option});
}

/**
 * Adjust left margin of chart
 * @param val
 */
function updateChartLeft(val) {
    document.getElementById("chartleftoutput").value = val + "%";
    var grid_option = chart.getOption().grid;
    _.forEach(grid_option, function (element) {
        element.left = val + "%";
    });
    chart.setOption({grid: grid_option});
}

/**
 * Adjust right margin of chart
 * @param val
 */
function updateChartRight(val) {
    document.getElementById("chartrightoutput").value = val + "%";
    var grid_option = chart.getOption().grid;
    _.forEach(grid_option, function (element) {
        element.right = val + "%";
    });
    chart.setOption({grid: grid_option});
}

/**
 * Adjust top margin of chart
 * @param val
 */
function updateChartTop(val) {
    document.getElementById("charttopoutput").value = val + "%";
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
 * Toggle to Amplicon Depth Label
 */
$("#toggleamplicondepthlabel").change(function () {
    var series_option = chart.getOption().series;
    if ($(this).prop("checked")) {
        _.forEach(series_option, function (element) {
            if (element.type === 'custom') {
                element.label.show = true
            }
        })
    } else {
        _.forEach(series_option, function (element) {
            if (element.type === 'custom') {
                element.label.show = false
            }
        })
    }
    chart.setOption({series: [...series_option]});
});

/**
 * Toggle to show gene label or not
 */
$("#togglegenelabel").change(function () {
    var series_option = chart.getOption().series;
    if ($(this).prop("checked")) {
        invisible_gene_label = false
    } else {
        invisible_gene_label = true
    }
    series_option[series_option.length - 1]["renderItem"] = renderGeneFeatures; // Re-update Gene Feature Chart Only
    chart.setOption({series: [...series_option]});
});

/**
 * Adjust the height of gene feature charts
 * @param val
 */
function updateGeneFeatureHeight(val) {
    document.getElementById("genefeatureheightoutput").value = val + "%";
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
                element.max = 40000;
            }
        });
        chart.setOption({yAxis: yaxis_option});
        document.getElementById("ymax").value = 40000;
    } else {
        _.forEach(yaxis_option, function (element) {
            if (element.gridIndex < yaxis_option.length - 1) {
                element.type = scale_type;
                element.min = 1;
                element.max = 100000;
            }
        });
        chart.setOption({yAxis: yaxis_option});
        document.getElementById("ymax").value = 100000;
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
 * Toggle tooltip for coverage chart, this action does not affect gene feature chart
 */
$("#toggletooltip").change(function () {
    if ($(this).prop("checked")) {
        chart.setOption({
            tooltip: {showContent: true},
        });
    } else {
        chart.setOption({
            tooltip: {showContent: false},
        });
    }
});

/**
 * Toggle slider zoom
 */
$("#toggleslider").change(function () {
    if ($(this).prop("checked")) {
        chart.setOption({
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
        });
    } else {
        chart.setOption({
            dataZoom: [
                {
                    type: "inside",
                    filterMode: "none",
                    xAxisIndex: [...Array(grid_length + 1).keys()],
                },
                {
                    show: false,
                    type: "slider",
                },
            ],
        });
    }
});

/**
 * Set range of gene feature to zoom in
 */
function setDataZoom() {
    var start = document.getElementById("start_pos").value;
    var end = document.getElementById("end_pos").value;
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
    document.getElementById("start_pos").value = 1;
    document.getElementById("end_pos").value = ref_len;
    chart.dispatchAction({
        type: "dataZoom",
        start: 0,
        end: 100,
    });
}

/**
 * Click gene feature to zoom in the view for this gene feature
 */
chart.on("click", function (params) {
    document.getElementById("start_pos").value = params.value[1];
    document.getElementById("end_pos").value = params.value[2];
    chart.dispatchAction({
        type: "dataZoom",
        startValue: params.value[1],
        endValue: params.value[2],
    });
});

/**
 * Double click to reset zoom to default
 */
chart.on("dblclick", function (params) {
    document.getElementById("start_pos").value = 1;
    document.getElementById("end_pos").value = ref_len;
    chart.dispatchAction({
        type: "dataZoom",
        start: 0,
        end: 100,
    });
});

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
    document.getElementById("chartheightinput").value = _.toInteger(height1);
    document.getElementById("chartheightoutput").value = _.toInteger(height1) + "%";
    document.getElementById("chartleftinput").value = _.toInteger(left);
    document.getElementById("chartleftoutput").value = _.toInteger(left) + "%";
    document.getElementById("chartrightinput").value = _.toInteger(right);
    document.getElementById("chartrightoutput").value = _.toInteger(right) + "%";
    document.getElementById("charttopinput").value = _.toInteger(top);
    document.getElementById("charttopoutput").value = _.toInteger(top) + "%";
    document.getElementById("genefeatureheightinput").value = _.toInteger(height2);
    document.getElementById("genefeatureheightoutput").value = _.toInteger(height2) + "%";
    // Set Axis to Log scale
    document.getElementById("scale").value = "log"
    document.getElementById("ymax").value = 100000
    // Reset Data Zoom
    resetDataZoom()
}

/**
 * Dispatch click, double click action on gene feature chart when render env changed or dark mode toggled
 */
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

/**
 * Write tooltip information to HTML table
 * @param headers, rows, classes
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

window.toTableHtml = toTableHtml;

/**
 * Calculate mean coverage for specific range
 * @param depths, start, end, gridIndex
 * @returns {number}
 */
function meanCoverage(depths, start, end, gridIndex) {
    var sub_array = _.slice(depths[gridIndex], start - 1, end);
    return _.mean(sub_array);
}

/**
 *  Calculate genome coverage depth which is >= low (10)
 * @param depths, start, end, gridIndex, low
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
 * @param depths, start, end, gridIndex
 * @returns {number}
 */
function medianCoverage(depths, start, end, gridIndex) {
    var sub_array = _.slice(depths[gridIndex], start - 1, end);
    return median(sub_array);
}