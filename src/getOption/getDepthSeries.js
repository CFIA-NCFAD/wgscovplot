import {getCoordsInterval} from "../util";

/**
 *  Plot mark area for low coverage regions
 *  @param {string} sample - Sample name
 *  @param {Array<string>} segments - An array of segments names
 *  @param {Array<Array<number>>} segmentsInterval - An array of segment start, end
 *  @param {} depths - Depths array for non segment virus or object for segment virus
 *  @param {number} lowCoverageThreshold - Low coverage threshold
 *  @param {boolean} nonVariantSites - whether to show non-variant sites information (tooltips options)
 *  @returns {Object} - Options for Mark Area
 *
 * */
function getMarkArea(sample, segments, segmentsInterval, depths, lowCoverageThreshold, opacity = 0.4) {
    let data = [];
    if (segments.length > 0) {
        for (let i = 0; i < segments.length; i++) {
            if (depths[sample][segments[i]] !== undefined && depths[sample][segments[i]] !== null) {
                let coords = getCoordsInterval(depths[sample][segments[i]], lowCoverageThreshold);
                for (let j = 0; j < coords.length; j++) {
                    let start = coords[j][0];
                    let end = coords[j][1];
                    data.push([
                        {
                            name: `Pos:${start} - ${end} (< ${lowCoverageThreshold}X)`,
                            xAxis: parseInt(start) + segmentsInterval[i][0] - 1,
                        },
                        {
                            xAxis: parseInt(end) + segmentsInterval[i][0] - 1,
                        }
                    ]);
                }
            }
        }
    } else {
        let coords = getCoordsInterval(depths, lowCoverageThreshold);
        for (let j = 0; j < coords.length; j++) {
            let start = coords[j][0];
            let end = coords[j][1];
            data.push([
                {
                    name: `Pos:${start} - ${end} (< ${lowCoverageThreshold}X)`,
                    xAxis: parseInt(start),
                },
                {
                    xAxis: parseInt(end),
                }
            ]);
        }
    }
    return {
        itemStyle: {
            color: "yellow",
            opacity: opacity
        },
        label: {
            show: false,
            position: "insideTop",
            fontSize: 10,
            rotate: 30,
            overflow: "truncate",
            ellipsis: "..."
        },
        data: data
    };
}

/**
 *  Plot mark lines for low coverage threshold
 *  @param {Array<Array<number>>} segmentsInterval - An array of segment start, end
 * @param {number} lowCoverageThreshold - Low coverage threshold
 * */
function getCoverageThresholdLine(segmentsInterval, lowCoverageThreshold) {
    return {
        silent: true,
        symbol: ["none", "none"],
        label: {
            show: true,
            formatter: '{c}'+'X',
        },
        lineStyle: {
            color: "#000",
            width: 1,
            type: "dashed",
            opacity: 0.4
        },
        data: [
            {
                name: "Low Coverage Threshold",
                yAxis: lowCoverageThreshold
            }
        ]
    };
}

/**
 * Define options for depth coverage charts
 * @param {Array<string>} samples - An array of samples names
 * @param {Array<string>} segments - An array of segments names
 * @param {} depths - Depths array for non segment virus or object for segment virus
 * @param {Array<Array<number>>} segmentsInterval - An array of segment start, end
 * @param {number} lowCoverageThreshold - Low coverage threshold
 * @param {boolean} nonVariantSites - whether to show tooltips for non-variant sites
 * @returns {Array<Object>}
 */
function getDepthSeries(samples, segments, depths, lowCoverageThreshold, segmentsInterval, nonVariantSites) {
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
            markLine: getCoverageThresholdLine(segmentsInterval, lowCoverageThreshold),
            markArea: (segments.length > 0) ? getMarkArea(samples[i], segments, segmentsInterval, depths, lowCoverageThreshold) :
                getMarkArea(samples[i], [], [], depths[i], lowCoverageThreshold),
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

export {getDepthSeries, getMarkArea, getCoverageThresholdLine};