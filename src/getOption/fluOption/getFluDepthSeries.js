/**
 *  Plot mark area for low coverage regions
 *  @param {string} sample - Sample name
 *  @param {Array<string>} segments - An array of segments names
 *  @param {Array<Array<number>>} segmentsRange - An array of segment start, end
 *  @param {Object} lowCoverageRegion - The object of low coverage regions
 *  @param {boolean} nonVariantSites - whether to show non-variant sites information (tooltips options)
 *  @returns {Object} - Options for Mark Area
 *
 * */
function getMarkArea(sample, segments, segmentsRange, lowCoverageRegion) {
    let data = [];
    for (let i = 0; i < segments.length; i++) {
        if (lowCoverageRegion[sample][segments[i]] !== '') {
            let coords = lowCoverageRegion[sample][segments[i]].split("; ");
            for (let j = 0; j < coords.length; j++) {
                let coord = coords[j].split("-");
                let start;
                let end;
                if (coord.length > 1) {
                    start = coord[0];
                    end = coord[1];
                } else { // single position
                    start = coord[0];
                    end = coord[0];
                }
                data.push([
                    {
                        name: `Region: ${start} - ${end}`,
                        xAxis: parseInt(start) + segmentsRange[i][0] - 1
                    },
                    {
                        xAxis: parseInt(end) + segmentsRange[i][0] - 1
                    }
                ]);
            }
        }
    }
    return {
        itemStyle: {
            color: 'yellow',
            opacity: 0.4
        },
        label: {
            show: false,
            position: 'insideTop',
            fontSize: 10,
            rotate: 30,
            overflow: 'truncate',
            ellipsis: '...'
        },
        data: data
    };
}

/**
 *  Plot mark lines to separate segments
 *  @param {Array<Array<number>>} segmentsRange - An array of segment start, end
 *
 * */
function getMarkLine(segmentsRange) {
    let data = [];
    for (let i = 0; i < segmentsRange.length; i++) {
        if (i === 0) {
            data.push({xAxis: segmentsRange[i][1]});
        } else if (i === segmentsRange.length - 1) {
            data.push({xAxis: segmentsRange[i][0]});
        } else {
            data.push({xAxis: segmentsRange[i][0]});
            data.push({xAxis: segmentsRange[i][1]});
        }
    }
    return {
        silent: true,
        symbol: ['none', 'none'],
        label: {
            show: false,
        },
        lineStyle: {
            color: "#000",
            width: 1,
            type: 'dashed',
            opacity: 0.5
        },
        data: data
    };
}

/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {Object} lowCoverageRegion - The object of low coverage regions
 * @param {Array<Array<number>>} segmentsRange - An array of segment start, end
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
            markLine: getMarkLine(segmentsRange),
            markArea: getMarkArea(samples[i], segments, segmentsRange, lowCoverageRegion),
            lineStyle: {
                color: "#666",
                opacity: 0,
            },
            tooltip: {
                trigger: nonVariantSites ? "axis" : "none"
            },
            silent: true,
            large: true,
        });
    }
    return depthSeries;
}

export {getFluDepthSeries};