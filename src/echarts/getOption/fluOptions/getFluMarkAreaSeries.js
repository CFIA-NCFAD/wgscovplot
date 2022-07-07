/**
 * Define options for mark area
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} depths - Object of depths array
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getFluMarkAreaSeries(samples, segments, depths) {
    let depthSeries = [];
    for (let i = 0; i < samples.length; i++) {
        for (let j =0; j < segments.length; j++){
            depthSeries.push({
                type: "line",
                xAxisIndex: i,
                yAxisIndex: i,
                markArea: {
                    itemStyle:{
                        color: 'rgba(255, 173, 177, 0.4)'
                    }
                },
                data:[
                    {
                        name: 'Segment',
                        xAxis:1
                    },
                    {
                        xAxis:1000
                    }
                ],
                silent: true,
            });
        }
    }
    return depthSeries;
}

export {getFluMarkAreaSeries}