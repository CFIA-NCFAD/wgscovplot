function getMarkArea(sample, segments, segmentsRange, lowNoCoverageRegions){
    let data = []
    for (let i = 0; i < segments.length; i++){
        if (lowNoCoverageRegions[sample][segments[i]] != ''){
            let coords= lowNoCoverageRegions[sample][segments[i]].split("; ")
            for (let j = 0; j < coords.length; j++){
                let coord = coords[j].split("-")
                let start;
                let end;
                if (coord.length > 1){
                    start = coord[0]
                    end = coord[1]
                }
                else { // single position
                    start = coord[0]
                    end = coord[0]
                }
                data.push([
                    {
                        name: `Region: ${start} - ${end}`,
                        xAxis: parseInt(start)+ segmentsRange[i][0] - 1
                    },
                    {
                        xAxis: parseInt(end) + segmentsRange[i][0] - 1
                    }
                ])
            }
        }
    }
    return {
        itemStyle:{
            color: 'yellow',
            opacity: 0.4
        },
        label: {
            show:false,
            position: 'insideTop',
            fontSize: 10,
            overflow: 'truncate',
            ellipsis:'...'
        },
        data: data
    }
}

function getMarkLine(segmentsRange){
    let data = []
    for (let i = 0; i < segmentsRange.length; i++){
        data.push({xAxis: segmentsRange[i][0]})
        data.push({xAxis: segmentsRange[i][1]})
    }
    return {
        symbol: ['none', 'none'],
        label: { show: false },
        data: data
    }
}

/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} lowCoverageRegion - The object of low coverage regions
 * @param {boolean} nonVariantSites - whether to show non-variant sites information (tooltips options)
 * @returns {Object} - Coverage Chart Option
 *
 * depths: { 'SAMPLE_NAME':{
 *                  'SEGMENT_NAME': []
 *              }
 *          }
 */
function getFluDepthSeries(samples, segments, lowCoverageRegion, segmentsRange, nonVariantSites) {
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
            //markLine: getMarkLine(segmentsRange),
            markArea: getMarkArea(samples[i], segments, segmentsRange, lowCoverageRegion),
            lineStyle: {
                color: "#666",
                opacity: 0,
            },
            tooltip:{
                trigger: nonVariantSites ? "axis" : "none"
            },
            silent: true,
            large: true,
        });
    }
    return depthSeries;
}

export {getFluDepthSeries}