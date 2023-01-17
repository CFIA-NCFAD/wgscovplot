import {get, isNil} from "lodash";
import {getCoordsInterval, NT_COLOURS} from "../util";
import {graphic} from "echarts/core";

/**
 *  Get options config for low coverage region highlighting using ECharts MarkArea
 *  @param {string} sample - Sample name
 *  @param {WgsCovPlotDB} db
 *  @returns {Object} - Options for Mark Area
 *
 * */
function getMarkArea({sample, db}) {
    const {
        depths,
        low_coverage_threshold = 10,
        lowCovAreaOpacity = 0.4,
        selectedSegments = null,
        segCoords = null,
    } = db;
    let data = [];
    if (!isNil(selectedSegments) && !isNil(segCoords) && selectedSegments.length > 0) {
        for (let i = 0; i < selectedSegments.length; i++) {
            let segment = selectedSegments[i];
            let sampleSegDepths = depths[sample][segment];
            if (isNil(sampleSegDepths)) {
                continue;
            }
            for (let {start, end} of getCoordsInterval(sampleSegDepths, low_coverage_threshold)) {
                data.push([
                    {
                        name: `${start}-${end} (<${low_coverage_threshold}X)`,
                        xAxis: start + segCoords[segment].start - 1,
                    },
                    {
                        xAxis: end + segCoords[segment].start - 1,
                    }
                ]);
            }
        }
    } else {
        for (let {start, end} of getCoordsInterval(depths[sample], low_coverage_threshold)) {
            data.push([
                {
                    name: `${start}-${end} (<${low_coverage_threshold}X)`,
                    xAxis: start,
                },
                {
                    xAxis: end,
                }
            ]);
        }
    }
    return {
        itemStyle: {
            color: "yellow",
            opacity: lowCovAreaOpacity
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
 * @param {WgsCovPlotDB} db
 * */
function getCoverageThresholdLine(db) {
    return {
        silent: true,
        symbol: ["none", "none"],
        label: {
            show: true,
            formatter: "{c}X",
        },
        lineStyle: {
            color: "#6b6464",
            width: 1,
            type: "dashed",
            opacity: db.showLowCovRegionsOpacity
        },
        data: [
            {
                name: "Low Coverage Threshold",
                yAxis: db.low_coverage_threshold
            }
        ]
    };
}

/**
 * Define options for depth coverage charts
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @returns {Object[]}
 */
function getDepthSeries(db) {
    const {
        selectedSamples,
        showNonVariantSiteTooltips,
    } = db;
    let depthSeries = [];
    for (let i = 0; i < selectedSamples.length; i++) {
        let sample = selectedSamples[i];
        depthSeries.push({
            markLine: getCoverageThresholdLine(db),
            markArea: getMarkArea({sample,db}),
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
            tooltip: {
                trigger: showNonVariantSiteTooltips ? "axis" : "none",
            },
            silent: true,
            large: true,
        });
    }
    return depthSeries;
}

/**
 * Define options for variant bar charts
 * @param {WgsCovPlotDB} db
 * @returns {Object[]}
 */
function getVariantsSeries(db) {
    const {
        selectedSamples,
        variants,
        depths,
        ref_seq,
        showVariantLabels,
        showVariantSiteTooltips = true,
        hideOverlappingVariantLabels = true,
    } = db;
    let variantSeries = [];
    let i = 0;
    for (let sample of selectedSamples) {
        let sampleVariants = get(variants, sample, []);
        variantSeries.push({
            type: "bar",
            xAxisIndex: i,
            yAxisIndex: i,
            data: sampleVariants.map((x) => [
                parseInt(x.POS),
                depths[sample][parseInt(x.POS) - 1],
            ]),
            barWidth: 2,
            itemStyle: {
                color: function ({data: [pos]}) {
                    let nt = ref_seq[pos - 1];
                    if (Object.prototype.hasOwnProperty.call(NT_COLOURS, nt)) {
                        return NT_COLOURS[nt];
                    }
                    return "#333";
                },
            },
            label: {
                show: showVariantLabels,
                position: "bottom",
                align: "left",
                verticalAlign: "middle",
                distance: 10,
                color: "inherit",
                rotate: -30,
                formatter: function ({data: [pos]}) {
                    let output = "";
                    Object.values(sampleVariants).forEach(({POS, REF, ALT}) => {
                        if (parseInt(POS) === pos) {
                            output += `${REF}${POS}${ALT}`;
                        }
                    });
                    return output;
                }
            },
            labelLayout: {
                hideOverlap: hideOverlappingVariantLabels
            },
            tooltip: {
                trigger: showVariantSiteTooltips ? "axis" : "none"
            }
        });
        i++;
    }
    return variantSeries;
}


/**
 * Renders the amplicon features shape
 * @param {boolean} show_amplicons - whether to plot amplicon feature or not (true for false)
 * @returns The shape of gene feature and/or amplicon feature
 */
function getRegionAmpliconDepthRenderer(show_amplicons) {
    /**
     *
     * @param params - Echarts params
     * @param api - Echarts api
     */
    function renderRegionAmpliconDepth({coordSys}, api) {
        let [startX, startY] = api.coord([api.value(0), api.value(2)]);
        let [endX, endY] = api.coord([api.value(1), 1]);
        let rectShape = graphic.clipRectByRect(
            {
                x: startX,
                y: startY,
                width: endX - startX,
                height: endY - startY
            },
            coordSys,
        );
        return rectShape && {
            type: "rect",
            shape: rectShape,
            style: api.style(),
            invisible: !show_amplicons
        };
    }

    return renderRegionAmpliconDepth;
}

/**
 * Define options for amplicon depth coverage bars
 * @param {WgsCovPlotDB} db
 * @returns {Object[]}
 */
function getRegionAmpliconDepthSeries(db) {
    const {
        selectedSamples,
        amplicon_depths,
        show_amplicons,
    } = db;
    let ampliconDepthSeries = [];
    for (let [i, sample] of selectedSamples.entries()) {
        ampliconDepthSeries.push({
            type: "custom",
            xAxisIndex: i,
            yAxisIndex: i,
            renderItem: getRegionAmpliconDepthRenderer(show_amplicons),
            label: {
                show: false,
                position: "top",
                distance: 25,
                rotate: 60
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
            data: amplicon_depths[sample],
        });
    }
    return ampliconDepthSeries;
}

export {
    getMarkArea,
    getCoverageThresholdLine,
    getDepthSeries,
    getVariantsSeries,
    getRegionAmpliconDepthSeries,
    getRegionAmpliconDepthRenderer,
};