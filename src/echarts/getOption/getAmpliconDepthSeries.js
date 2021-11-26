import {graphic} from "echarts/core";

/**
 * Define options for amplicon depth coverage bars
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Dict[]>} ampliconDepthBarData - Array of dictionary geneFeature or amplicon data
 * @returns {Array<Dict[]>}
 */
function getAmpliconDepthSeries(samples, ampliconDepthBarData) {
    var ampliconDepthSeries = [];
    if (amplicon === 'True'){
        for (var [i, sample] of samples.entries()) {
            ampliconDepthSeries.push({
                type: "custom",
                xAxisIndex: i,
                yAxisIndex: i,
                renderItem: function (params, api) {
                    var start = api.coord([api.value(0), api.value(2)]);
                    var end = api.coord([api.value(1), 1]);
                    var rectShape = graphic.clipRectByRect(
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
                data: ampliconDepthBarData[sample],
            })
        }
    }
    return ampliconDepthSeries;
}

export {getAmpliconDepthSeries};