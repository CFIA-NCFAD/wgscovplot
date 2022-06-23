/**
 * Initialize the events handler for chart
 */
function initWgscovplotEvent(){
    $(document).ready(function () {


        $("#selectedsamples").select2({
            tags: true,
        });

        /**
        * Jquery actions to make the list of samples is not forced in alphabetical order
        */
        $("#selectedsamples").on("select2:select", function (evt) {
            let element = evt.params.data.element;
            let $element = $(element);
            $element.detach();
            $(this).append($element);
            $(this).trigger("change");
        });

        $("#selectedsegments").select2({
            tags: true,
        });

        $("#selectedsamples").on("change", function () {
            updateFluCoverageChartOption(chart.getOption());
        });

        $("#selectedsegments").on("change", function () {
            updateFluCoverageChartOption(chart.getOption());
        });

    });
}


function updateFluCoverageChartOption(charOption){
    const [plotSamples, plotSegments] = getCurrentSamplesSegments(charOption)
    let updateOption = wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, window.depths, window.variants, window.ref_seq)
    chart.setOption(option=updateOption, {notMerge:true})
    chart.resize({
        width: 'auto'
    })
}

/**
 * When the chart is initialized, the first 3 samples, segments are plotted
 * @param {Array<string>} samples - An array of samples name
 */
function setDefaultSamplesSegments(samples, segments) {
    // Set default samples display
    let $selectedsamples = $("#selectedsamples");
    $selectedsamples.select2();
    $selectedsamples.val(samples);
    $selectedsamples.trigger('change');

    // Set default samples display
    let $selectedsegments= $("#selectedsegments");
    $selectedsegments.select2();
    $selectedsegments.val(segments);
    $selectedsegments.trigger('change');
}

/**
 * Get the list of current samples on the control menu
 * @returns {*[][]} An array of samples name
 */
function getCurrentSamplesSegments(chartOption) {
    let samples = [];
    let segments = [];
    // If the chart is not initialized yet, get 3 first samples from window.samples
    if (chartOption === undefined || chartOption === null){
        samples = window.samples.slice(0, 3);
        segments = window.segments.slice(0, 8);
    } else{
        let selectData1 = $("#selectedsamples").select2("data");
        let selectData2 = $("#selectedsegments").select2("data");
        for (let [key1, entries1] of selectData1.entries()) {
            samples.push(selectData1[key1].text);
        }
        for (let [key2, entries2] of selectData2.entries()) {
            segments.push(selectData2[key2].text);
        }
    }
    return [samples, segments];
}

function initWgscovplotRenderEnv(){

    let chartOption = chart.getOption();
    const [plotSamples, plotSegments] = getCurrentSamplesSegments(chartOption);
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamplesSegments(plotSamples, plotSegments);
        chart.setOption(option=wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, window.depths, window.variants, window.ref_seq))
    }
    chart.setOption(option=option)
    onChartDataZoomActions()
}

/**
 * Set scale for y Axis
 */
function setScale() {
    let scaleType = document.getElementById("scale").value;
    let yAxisOption = chart.getOption().yAxis;
    yAxisOption.forEach(element => {
        if (element.gridIndex < yAxisOption.length -1) {
            element.type = scaleType;
        }
    });
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
 * Dispatch click/dbclick actions for the whole chart
 */
function onChartDataZoomActions(){

    chart.on("click", function(params){
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === 'custom'){
            setDataZoom(params.value.start, params.value.end);
        }

    });

    chart.on("dblclick", function(params){
        if (params.componentIndex === chart.getOption().series.length - 1 && params.componentSubType === 'custom'){
            let chartOption = chart.getOption();
            const [plotSamples, plotSegments] = getCurrentSamplesSegments(chartOption);
            const maxSegmentsLength = wgscovplot.getMaxSegmentsLength(plotSamples, plotSegments, window.depths)
            const xAxisMax = wgscovplot.getXAxisMax(maxSegmentsLength);
            setDataZoom(1, xAxisMax);
        }
    });

}

/**
 * Set zoom view for the chart
 * @param {number} zoomStart - Start view point
 * @param {number} zoomEnd - End view point
 */
function setDataZoom(zoomStart, zoomEnd){
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
                (parseFloat(gridOption[i - 1].top.replace("%","")) +
                parseFloat(gridOption[i - 1].height) +
                parseFloat(val)).toFixed(1) + "%";
        }
    }
    chart.setOption({grid: gridOption});
}

/**
 * Reset Grid Dislay to optimal configuration
 */
function resetGridDisplay(){
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
        document.getElementById("chart-left-output").value = left+ "%";

        document.getElementById("chart-right-input").value = right;
        document.getElementById("chart-right-output").value = right + "%";
    }
}
