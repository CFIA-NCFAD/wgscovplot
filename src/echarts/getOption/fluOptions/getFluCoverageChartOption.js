import {getFluDatasets} from "./getFluDataSets";
import {getFluYAxes} from "./getFluAxes";
import {getFluXAxes} from "./getFluAxes";
import {getFluGrids} from "./getFluGrids";
import {getFluDepthSeries} from "./getFluDepthSeries";

function getFluCoverageChartOption(samples, segments, depths) {
    let chartOptions = {
        title: {},
        dataset: getFluDatasets(depths),
        xAxis: getFluXAxes(samples, segments, depths),
        yAxis: getFluYAxes(samples, segments),
        series: [
            ...getFluDepthSeries(samples, segments)
        ],
        tooltip:{
            trigger: "axis",
            enterable: true,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            position: "cursor",
            axisPointer: {
                type: 'line'
            },
        },
        grid: getFluGrids(samples, segments)
    };
    return chartOptions;
}

export {getFluCoverageChartOption};