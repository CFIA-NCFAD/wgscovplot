/**
 * Get depths and variants for samples
 * @param {Array<string>} samples - An array of samples name
 * @param {Object} depthsData - Depths data
 * @param {Object} variantsData - Variants data
 * @returns [Array<Array<number>>, Array<Array<Object>>]
 */
function getDepthsVariants(samples, depthsData, variantsData){
    let depths = [];
    let variants = [];
    samples.forEach(sample => {
        if (depthsData[sample] !== undefined && depthsData[sample] !== null) {
            depths.push(depthsData[sample]);
        } else{
            depths.push([]);
        }
        if (variantsData[sample] !== undefined && variantsData[sample] !== null){
            variants.push(variantsData[sample]);
        }else{
            variants.push({});
        }
    });
    return [depths, variants];
}
/**
 * Updates options for coverage charts.
 * Whenever the selected samples changed, chart options such as y Axis scale, yMax, DataZoom are reserved
 * Users's settings are respected by keeping old settings and set it back.
 * @param {Array<string>} samples - An array of samples name
 */
function updateCoverageChartOption(samples) {
    let scaleType;
    let yAxisMax;
    const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
    let chartOption = chart.getOption();
    // Reserver Tooltip option
    let isTooltipEnable = document.getElementById("toggle-tooltip").checked;
    let triggerOnType;
    if (isTooltipEnable){
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    }else{
        triggerOnType ="none";
    }
    let isVariantSites = document.getElementById("toggle-tooltip-variant-sites").checked;
    let isNonVarianSites = document.getElementById("toggle-tooltip-non-variant-sites").checked;
    let isVariantComparison = document.getElementById("toggle-variant-comparison").checked;
    let isCoverageStatView = document.getElementById("toggle-coverage-stat").checked;
    let isShowMutation = document.getElementById("toggle-mutation").checked;
    let isShowXAxisLabel = document.getElementById("toggle-xaxis-label").checked;
    let isHideOverlapMutation = document.getElementById("toggle-hideoverlap-mutation").checked;
    let updateOption = wgscovplot.getCoverageChartOption(geneAmpliconFeatureData, regionAmpliconDepthData, window.refSeq,
        yAxisMax, samples, depths, variants, geneFeature, amplicon,
        triggerOnType= triggerOnType, isVariantSites=isVariantSites,
        isNonVariantSites=isNonVarianSites, isInfoComparison=isVariantComparison,
        isCovergateStatView=isCoverageStatView,
        isShowMutation=isShowMutation,
        isShowXAxisLabel=isShowXAxisLabel,
        isHideOverlapMutation=isHideOverlapMutation);

    // Reserve tooltip in series option
    let seriesOption = updateOption.series;
    seriesOption.forEach(element => {
        if (element.type === 'line'){
            element.tooltip.trigger = isNonVariantSites ? "axis" : "none";
        }
        else if (element.type === 'bar'){
            element.tooltip.trigger = isVariantSites ? "axis" : "none";
        }
    });
    let isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
    updateOption.tooltip[0].position = tooltipPosition(isFixTooltipPostion);

    // Reserve datzoom
    let oldDataZoom = chartOption.dataZoom;
    updateOption.dataZoom = oldDataZoom;
    updateOption.dataZoom.forEach(element => {
        element.xAxisIndex = [...Array(updateOption.grid.length).keys()];
    });

    //Reserve scale and yAxis max
    scaleType = $("#scale").val();
    yAxisMax = $("#ymax").val();
    updateOption.yAxis = updateYAxisOption(updateOption.yAxis, scaleType, yAxisMax);

    //Toggle amplicon
    let isShowAmplicon;
    if (amplicon){
        isShowAmplicon = document.getElementById("toggle-amplicon").checked;
        for (var i = 0; i < seriesOption.length - 1; i++){
            if (seriesOption[i].type === "custom"){
                seriesOption[i].renderItem = wgscovplot.getRegionAmpliconDepthRenderer(isShowAmplicon);
            }
        }
        updateOption.grid = updateGrid(geneFeature, isShowAmplicon);
    }

    if (geneFeature){
        isShowAmplicon = amplicon ? document.getElementById("toggle-amplicon").checked : amplicon;
        seriesOption[seriesOption.length - 1].renderItem = wgscovplot.getGeneFeatureRenderer(document.getElementById("toggle-genelabel").checked,
            geneAmpliconFeatureData, isShowAmplicon);
    }

        // Reserve grid option
    updateOption.grid.forEach(element => {
        element.left = $("#chart-left-input").val() + "%";
        element.right = $("#chart-right-input").val() + "%";
    });

    //set chart option
    chart.setOption(option=updateOption, {notMerge:true});
    // Update control menu
    updateControlMenu();
    variantHeatmap.setOption(option=wgscovplot.getVariantHeatmapOption(samples, window.variants));
}

/**
 * Update scale type and max for Y Axis
 * @param {Object} yAxisOption - Options of Yaxis need to be updated
 * @param {string} scaleType - Either log or value
 * @param {number} yAxisMax - Max value is set for Y Axis
 * @returns {Object} Returns the updated options (Scale type or ymax) for yAxis
 */
function updateYAxisOption(yAxisOption, scaleType, yAxisMax){
    let len = (amplicon || geneFeature) ? yAxisOption.length - 1 : yAxisOption.length;
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
    $selectedsamples.trigger('change');
}

/**
 * Get the list of current samples on the control menu
 * @returns {Array<string>} An array of samples name
 */
function getCurrentSamples(chartOption) {
    let samples = [];
    // If the chart is not initalizaed yet, get 3 first samples from window.samples
    if (chartOption === undefined || chartOption === null){
        samples = window.samples.slice(0, 3);
    } else{
        let selectData = $("#selectedsamples").select2("data");
        for (let [key, entries] of selectData.entries()) {
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

        $("#selected-gene-feature-name").select2({
            tags: true,
            width: '100%',
        });


        /**
         * Update chart options when adding/romoving samples
         * The chart options such as y Axis scale, yMax, DataZoom are reserved (users's settings are respected)
         */
        $("#selectedsamples").on("change", function () {
            updateCoverageChartOption(getCurrentSamples(chart.getOption()));
        });

        $("#selected-gene-feature-name").on("change", function () {
            $("#selected-gene-feature-name").select2('data');
        });

        /**
         * Toggle to show gene label or not
         */
        $("#toggle-genelabel").change(function () {
            let seriesOption = chart.getOption().series;
            let showGeneLabel = $(this).prop("checked");
            let isShowAmplicon = amplicon ? document.getElementById("toggle-amplicon").checked : amplicon;
            seriesOption[seriesOption.length - 1].renderItem = wgscovplot.getGeneFeatureRenderer(showGeneLabel, geneAmpliconFeatureData, isShowAmplicon);
            chart.setOption({series: [...seriesOption]});
        });

        /**
         * Toggle to Amplicon Depth Plot
         */
        $("#toggle-amplicon").change(function () {
            let isChecked = $(this).prop("checked");
            let seriesOption = chart.getOption().series;
            for (var i = 0; i <seriesOption.length - 1; i++){
                if (seriesOption[i].type === "custom"){
                    seriesOption[i].renderItem = wgscovplot.getRegionAmpliconDepthRenderer(isChecked);
                }
            }
            if (geneFeature){
                let isShowGeneLabel = document.getElementById("toggle-genelabel").checked;
                seriesOption[seriesOption.length - 1].renderItem = wgscovplot.getGeneFeatureRenderer(isShowGeneLabel, geneAmpliconFeatureData, isChecked);
            }
            var gridOption = updateGrid(geneFeature, isChecked);
            gridOption.forEach(element => {
                element.left = $("#chart-left-input").val() + "%";
                element.right = $("#chart-right-input").val() + "%";
            });
            chart.setOption({series: [...seriesOption], grid: [...gridOption]});
            updateControlMenu();
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

        /**
         * Toggle to show X Axis
         */
        $("#toggle-xaxis-label").change(function () {
            let xAxisOption = chart.getOption().xAxis;
            let showAxisLabel = $(this).prop("checked");
            let gridLength = (amplicon || geneFeature) ? xAxisOption.length - 1 : xAxisOption.length;
            for (let i = 0; i < gridLength; i++){
                xAxisOption[i].axisLabel.show= showAxisLabel;
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
                tooltip: {triggerOn: isChecked ? triggerType: "none"},
            });
        });

                /**
         * Fix tooltip position on the charts
         * If it is fixed: tooltip will be on the right if mouse hovering on the left and vice versa
         * Default tooltip follows cursor
         */
        $("#fix-tooltip-postion").change(function () {
            let isChecked = $(this).prop("checked");
            let tooltipOption  = chart.getOption().tooltip;
            tooltipOption[0].position = tooltipPosition(isChecked);
            chart.setOption({tooltip: tooltipOption});
        });

        /**
         * Toogle tooltip trigger on click
         */
        $("#toggle-tooltip-trigger-click").change(function () {
            let isChecked = $(this).prop("checked");
            let isTooltipEnable = document.getElementById("toggle-tooltip").checked;
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
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                isChecked, document.getElementById("toggle-tooltip-non-variant-sites").checked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });

        $("#toggle-tooltip-non-variant-sites").change(function (){
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked, isChecked,
                document.getElementById("toggle-variant-comparison").checked,
                document.getElementById("toggle-coverage-stat").checked);
        });


        $("#toggle-variant-comparison").change(function (){
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let samples = getCurrentSamples(chartOption);
            const [depths, variants] = getDepthsVariants(samples, window.depths, window.variants);
            updateTooltipOption(samples, depths, variants, seriesOption,
                document.getElementById("toggle-tooltip-variant-sites").checked,
                document.getElementById("toggle-tooltip-non-variant-sites").checked,
                isChecked, document.getElementById("toggle-coverage-stat").checked);
        });

        $("#toggle-coverage-stat").change(function (){
            let isChecked = $(this).prop("checked");
            let chartOption = chart.getOption();
            let seriesOption = chartOption.series;
            let samples = getCurrentSamples(chartOption);
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
            let isChecked = $(this).prop("checked");
            let numChart = chart.getOption().grid.length;
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
            });
            variantHeatmap.setOption({
                dataZoom: [
            {
                type: "inside"
            },
            {
                type: "slider",
                show: isChecked
            }
        ]
            });
        });
    });
}

/**
 * Initialize the enviroment for the chart (sgv/canvas or dark/white mode)
 * The entire chart will be disposed and re-initialized
 * However, the old settings of charts are reserved (users's settings are respected)
 */
function initWgscovplotRenderEnv() {
    let chartOption = chart.getOption();
    let plotSamples = getCurrentSamples(chartOption);
    const [plotDepths, plotVariants] = getDepthsVariants(plotSamples, window.depths, window.variants);
    if (chartOption === undefined || chartOption === null) {
        setDefaultSamples(plotSamples);
        chart.setOption(option=wgscovplot.getCoverageChartOption(geneAmpliconFeatureData, regionAmpliconDepthData, window.refSeq,
            maxDepth, plotSamples, plotDepths, plotVariants, geneFeature, amplicon));
        variantHeatmap.setOption(option=wgscovplot.getVariantHeatmapOption(plotSamples, window.variants));
    } else {
        let renderEnv = document.getElementById("render-env").value;
        let isChecked = document.getElementById("toggle-darkmode").checked;
        let mode = isChecked ? "dark" : "white";
        let gridOption = chart.getOption().grid;
        let seriesOption = chart.getOption().series;
        let scaleType = chart.getOption().yAxis[0].type;
        let yAxisMax = chart.getOption().yAxis[0].max;
        let dataZoomOption = chart.getOption().dataZoom;
        let tooltipOption = chart.getOption().tooltip;
        wgscovplot.echarts.dispose(chart); // destroy chart instance and re-init chart
        $chart = document.getElementById("chart");
        chart = wgscovplot.echarts.init($chart, mode, {renderer: renderEnv});
        let option = wgscovplot.getCoverageChartOption(geneAmpliconFeatureData, regionAmpliconDepthData, window.refSeq,
            yAxisMax, plotSamples, plotDepths, plotVariants, geneFeature, amplicon);
        // Keep grid option
        option.grid = gridOption;
        // Keep data zoom option
        option.dataZoom = dataZoomOption;
        // Keep yAxis option
        option.yAxis = updateYAxisOption(option.yAxis, scaleType, yAxisMax);
        // Keep tooltip
        option.tooltip = tooltipOption;
        // Keep series
        option.series = seriesOption;
        //set chart option
        chart.setOption(option = option);
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
    let isTooltipEnable = document.getElementById("toggle-tooltip").checked;
    let triggerOnType;
    if (isTooltipEnable){
        triggerOnType = document.getElementById("toggle-tooltip-trigger-click").checked ? "click" : "mousemove";
    }else{
        triggerOnType ="none";
    }
    seriesOption.forEach(element => {
        if (element.type === 'line'){
            element.tooltip.trigger = isNonVariantSites ? "axis" : "none";
        }
        else if (element.type === 'bar'){
            element.tooltip.trigger = isVariantSites ? "axis" : "none";
        }
    });
    let tooltipOption = wgscovplot.getTooltips(samples, depths, variants, window.refSeq,
    triggerOnType=triggerOnType, isInfoComparison=isInfoComparison, isCovergateStatView=isCovergateStatView);
    let isFixTooltipPostion = document.getElementById("fix-tooltip-postion").checked;
    tooltipOption[0].position = tooltipPosition(isFixTooltipPostion);
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
            let obj = {top: 5};
            obj[['left', 'right'][+(pos[0] < size.viewSize[0] / 2)]] = 5;
            return obj;
        };
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
    let gridOption = variantHeatmap.getOption().grid;
    gridOption[0].height = val + "%";
    variantHeatmap.setOption({grid: gridOption});
}

/**
 * Apply View for selected gene feature
 */
function applyFeatureView(){
    let featureName = $("#selected-gene-feature-name").select2('data');
    let start = [], end = [];
    let minStart, maxEnd;
    featureName.forEach(x => {
        for (let i = 0; i < geneAmpliconFeatureData.length; i++){
            if (x.text === geneAmpliconFeatureData[i].name){
                start.push(geneAmpliconFeatureData[i].value.start);
                end.push(geneAmpliconFeatureData[i].value.end);
                break;
            }
        }
    });
    minStart = Math.min(...start);
    maxEnd = Math.max(...end);
    if (featureName.length === 0){
        minStart = 1;
        maxEnd = refSeqLength;
    }
    setDataZoom(minStart, maxEnd);
}

/**
 * Adjust subplot height
 * @param {number} val - Subplots height percent value
 */
function updateSubPlotHeight(val) {
    document.getElementById("chart-height-output").value = val + "%";
    let gridOption = chart.getOption().grid;
    let len = (amplicon || geneFeature) ? gridOption.length - 1 : gridOption.length;
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
    if (amplicon || geneFeature){
        gridOption[len].top =
                parseFloat(gridOption[len - 1].top) +
                parseFloat(gridOption[len - 1].height) +
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
 * Adjust the height of gene/amplicon feature charts
 * @param {number} val - Subplots gene feature height percent value
 */
function updateGeneFeatureHeight(val) {
    document.getElementById("genefeature-height-output").value = val + "%";
    let gridOption = chart.getOption().grid;
    gridOption[gridOption.length - 1].height = val + "%";
    chart.setOption({grid: gridOption});
}

/**
 * Set scale for y Axis
 */
function setScale() {
    let scaleType = document.getElementById("scale").value;
    let yAxisMax = document.getElementById("ymax").value;
    let yAxisOption =  updateYAxisOption(chart.getOption().yAxis, scaleType, yAxisMax)
    chart.setOption({yAxis: yAxisOption});
}

/**
 * Set yMax for Y Axis
 */
function setYMax() {
    let yMax = document.getElementById("ymax").value;
    let yAxisOption = chart.getOption().yAxis;
    let len = (amplicon || geneFeature) ? yAxisOption.length - 1 : yAxisOption.length;
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
    let start;
    let end;
    if (zoomStart === 0 && zoomEnd === 0){
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
 * Update grid configuration for chart
 * @param geneFeature
 * @param amplicon
 * @returns {Array<Object>} grid configuration for the charts
 */
function updateGrid(geneFeature, amplicon){
    let chartOption = chart.getOption();
    let currentSamples = getCurrentSamples(chartOption);
    let doubleStrand = false;
    for (let i = 0; i < geneAmpliconFeatureData.length; i++){
        if (geneAmpliconFeatureData[i].value.strand === -1){
            doubleStrand = true;
            break;
        }
    }
    let gridOption = wgscovplot.getGrids(currentSamples, geneFeature, amplicon, doubleStrand);
    return gridOption;
}

/**
 * Reset Grid Dislay to optimal configuration
 */
function resetGridDisplay(){
    let isShowAmplicon = amplicon ? document.getElementById("toggle-amplicon").checked : amplicon;
    let gridOption = updateGrid(geneFeature, isShowAmplicon);
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
        if (amplicon || geneFeature){
            document.getElementById("genefeature-height-input").value = parseFloat(gridOption[gridOption.length-1].height);
            document.getElementById("genefeature-height-output").value = parseFloat(gridOption[gridOption.length-1].height)+ "%";
        }
    }
}

