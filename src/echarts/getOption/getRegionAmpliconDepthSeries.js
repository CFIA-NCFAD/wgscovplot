import {graphic} from "echarts/core";

/**
 * Renders the amplicon features shape
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @returns The shape of gene feature and/or amplicon feature
 */
function getRegionAmpliconDepthRenderer (amplicon){
    /**
     *
     * @param params - Echarts params
     * @param api - Echarts api
     */
    function renderRegionAmpliconDepth(params, api){
        let start = api.coord([api.value(0), api.value(2)]);
        let end = api.coord([api.value(1), 1]);
        let rectShape = graphic.clipRectByRect(
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
            style: api.style(),
            invisible: !amplicon
        };
    }
    return renderRegionAmpliconDepth;
}
/**
 * Define options for amplicon depth coverage bars
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Object>} regionAmpliconDepthData - Array of dictionary geneFeature or amplicon data
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true or false)
 * @returns {Array<Object>}
 */
function getRegionAmpliconDepthSeries(samples, regionAmpliconDepthData, amplicon) {
    let ampliconDepthSeries = [];
    if (amplicon){
        for (let [i, sample] of samples.entries()) {
            ampliconDepthSeries.push({
                type: "custom",
                xAxisIndex: i,
                yAxisIndex: i,
                renderItem: getRegionAmpliconDepthRenderer(amplicon),
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
                tooltip: {
                    trigger: "none"
                },
                silent: true,
                data: regionAmpliconDepthData[sample],
            });
        }
    }
    return ampliconDepthSeries;
}

export {getRegionAmpliconDepthSeries, getRegionAmpliconDepthRenderer};