/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {boolean} nonVariantSites - whether to show non-variant sites information (tooltips options)
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getFluDepthSeries(samples, segments, nonVariantSites) {
    let depthSeries = [];
    for (let i = 0; i < samples.length; i++) {
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
            tooltip:{
                trigger: nonVariantSites ? "axis" : "none"
            },
            /*
            markArea: {
                itemStyle:{
                    color: 'rgba(255, 173, 177, 0.4)'
                },
                data:[
                    {
                        name: 'Segment',
                        xAxis:1
                    },
                    {
                        xAxis:1000
                    }
                ]
            },*/
            silent: true,
            large: true,
        });
    }
    return depthSeries;
}

export {getFluDepthSeries}