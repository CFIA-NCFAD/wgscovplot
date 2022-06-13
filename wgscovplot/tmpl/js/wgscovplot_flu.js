/**
 * Initialize the events handler for chart
 */
function initWgscovplotEvent(){
    $(document).ready(function () {

        /**
         * Jquery actions to make the list of samples is not forced in alphabetical order
         */
        $("#selectedsamples").select2({
            tags: true,
        });

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

        $("#selectedsegments").on("select2:select", function (evt) {
            let element = evt.params.data.element;
            let $element = $(element);
            $element.detach();
            $(this).append($element);
            $(this).trigger("change");
        });

        $("#selectedsamples").on("change", function () {
            updateFluCoverageChartOption(chart.getOption());
        });

        $("#selectedsegments").on("change", function () {
            updateFluCoverageChartOption(chart.getOption());
        });

    });
}

function getFluDepths (samples, segments){
    let depths = [];
    for (let i = 0; i < samples.length; i++){
        for (let j = 0; j < segments.length; j++){
            depths.push(window.depths[samples[i]][segments[j]])
        }
    }
    console.log("Number of depths", depths.length)
    return depths;
}

function updateFluCoverageChartOption(charOption){
    const [plotSamples, plotSegments] = getCurrentSamplesSegments(charOption)
    let plotDepths = getFluDepths(plotSamples, plotSegments);
    let updateOption = wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, plotDepths)
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
        segments = window.segments.slice(0, 3);
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
    let plotDepths = getFluDepths(plotSamples, plotSegments)
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamplesSegments(plotSamples, plotSegments);
        chart.setOption(option=wgscovplot.getFluCoverageChartOption(plotSamples, plotSegments, plotDepths))
    }
    chart.setOption(option=option)
}