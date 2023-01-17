
import {FEATURE_PLOT_PROPS, toTableHtml} from "./util";
import {graphic} from "echarts/core";

/**
 * Define the points of gene/amplicon features shape
 * @param {number} x - x-axis coordinate
 * @param {number} y - y-axis coordinate
 * @param {number} width - width of shape
 * @param {number} height - height of shape
 * @param {number} strand - strand of feature
 * @param {string} type - type of feature
 * @returns {(number[])[] | null} - Coordinate of 5 points for gene feature or 4 points for amplicon feature
 */
function shapePoints(
    {
        x,
        y,
        width,
        height,
        feature: {
            value: {
                strand,
                type,
            }
        }
    }) {
    let taperOffset;
    if (width < 10) { // Element width is too small
        taperOffset = width / 2;
    } else {
        taperOffset = 5;
    }
    if (type === "gene") {
        if (strand === 1) {
            return [
                [x, y],
                [x + width - taperOffset, y],
                [x + width, y - height / 2],
                [x + width - taperOffset, y - height],
                [x, y - height],
            ];
        } else {
            return [
                [x, y - height / 2],
                [x + taperOffset, y],
                [x + width, y],
                [x + width, y - height],
                [x + taperOffset, y - height],
            ];
        }
    } else if (type === "amplicon" || type === "segment") {
        return [
            [x, y],
            [x + width, y],
            [x + width, y - height],
            [x, y - height],
        ];
    } else {
        return null;
    }
}

/**
 * Renders the gene/amplicon features shape
 * @param {WgsCovPlotDB} db
 * @returns The shape of gene feature and/or amplicon feature
 */
function getGeneFeatureRenderer(db) {
    const {
        showGeneLabels,
        echart_features,
        show_amplicons,
    } = db;
    let yStart = 0;

    function renderGeneFeatures(params, api) {
        const {
            dataIndex,
            coordSys,
        } = params;
        const feature = echart_features[dataIndex];
        const {
            name,
            value: {
                start,
                end,
                level,
                rotate
            },
        } = feature;
        const leftCoord = coordSys.x;
        const rightCoord = coordSys.width + coordSys.x;
        const [startX, startY] = api.coord([start, dataIndex]);
        if (dataIndex === 0) {
            yStart = startY;
        }
        const [endX] = api.coord([end, dataIndex]);
        const height = FEATURE_PLOT_PROPS.rec_items_height;
        const width = endX - startX;
        const y = yStart - height / 2 - level;
        const points = shapePoints({
            x: startX,
            y,
            width,
            height,
            feature,
        });
        const shape = graphic.clipPointsByRect(points, coordSys);
        let invisible = false;
        if (feature.value.type === "gene") {
            if (showGeneLabels) {
                // Element width is too small and hide label at the edges
                invisible = width < 10 || startX >= rightCoord || endX <= leftCoord;
            } else {
                invisible = true;
            }

            return {
                type: "polygon",
                shape: {
                    points: shape,
                },
                style: api.style(),
                textContent: {
                    type: "text",
                    invisible: invisible,
                    style: {
                        text: name,
                        fill: feature.itemStyle.color,
                        fontStyle: "normal",
                        fontSize: 12,
                        fontWeight: "bolder",
                    },
                },
                textConfig: {
                    position: "top",
                    distance: 18,
                    rotation: rotate,
                    offset: "center",
                    local: true,
                },
            };
        } else if (feature.value.type === "amplicon") {
            return {
                type: "polygon",
                shape: {
                    points: shape,
                },
                style: api.style(),
                textContent: {},
                textConfig: {},
                invisible: !show_amplicons
            };
        } else if (feature.value.type === "segment") {
            invisible = width < 10 || startX >= rightCoord || endX <= leftCoord;
            return {
                type: "polygon",
                shape: {
                    points: shape,
                },
                style: api.style(),
                textContent: {
                    type: "text",
                    invisible: invisible,
                    style: {
                        text: name,
                        fill: feature.itemStyle.color,
                        fontStyle: "normal",
                        fontSize: 12,
                        fontWeight: "bolder",
                    },
                },
                textConfig: {
                    position: "top",
                    distance: 18,
                    rotation: rotate,
                    offset: "center",
                    local: true,
                },
            };
        } else {
            return null;
        }
    }

    return renderGeneFeatures;
}

/**
 * Define options for features charts
 * @param {WgsCovPlotDB} db - wgscovplot DB object
 * @param {number} index - gene feature is displayed in the last index of grid
 * @returns {Object}
 */
function getGeneFeatureSeries(
    {
        db,
        index,
    }) {
    const {
        echart_features,
    } = db;
    return {
        type: "custom",
        xAxisIndex: index,
        yAxisIndex: index,
        renderItem: getGeneFeatureRenderer(db),
        labelLayout: {
            hideOverlap: false,
        },
        data: echart_features,
        tooltip: {
            trigger: "item",
            enterable: true,
            appendToBody: true,
            triggerOn: "mousemove",
            renderMode: "html",
            borderRadius: 6,
            borderWidth: 2,
            showContent: "true",
            position: "top",
            textStyle: {
                fontSize: 15,
                fontWeight: "bolder",
            },
            formatter: function ({name, value: {start, end, strand}}) {
                let rows = [
                    ["Range", `${start.toLocaleString()}-${end.toLocaleString()}`],
                    ["Length", (end - start + 1).toLocaleString()],
                    ["Strand", strand],
                ];
                return toTableHtml({
                    headers: [name, ""],
                    rows,
                });
            },
        },
    };

}

export {getGeneFeatureSeries, getGeneFeatureRenderer};