import {shapePoints} from "./shapePoints";
import {geneFeaturePlotProperties} from "../util";
import {graphic} from "echarts/core";

/**
 * Renders the gene/amplicon features shape
 * @param {boolean} isShowGeneLabel - Whether to show gene label or not
 * @param {Array<Object>} geneFeatureAmpliconData - Array of dictionary geneFeature or amplicon data
 * @param {boolean} amplicon - whether to plot amplicon feature or not (true for false)
 * @returns The shape of gene feature and/or amplicon feature
 */
function getGeneFeatureRenderer(isShowGeneLabel, geneFeatureAmpliconData, amplicon) {
    /**
     *
     * @param params - Echarts params
     * @param api - Echarts api
     */
    let yStart;
    function renderGeneFeatures(params, api) {
        let points, shape, rotateAngle;
        let start, end, height, width, x, y, leftCoord, rightCoord;
        const categoryIndex = params.dataIndex;
        const feature = geneFeatureAmpliconData[categoryIndex];
        leftCoord = params.coordSys.x;
        rightCoord = params.coordSys.width + params.coordSys.x;
        start = api.coord([feature.value.start, categoryIndex]);
        if (categoryIndex === 0) {
            yStart = start[1];
        }
        end = api.coord([feature.value.end, categoryIndex]);
        height = geneFeaturePlotProperties.rec_items_height;
        width = end[0] - start[0];
        x = start[0];
        y = yStart - height / 2 - feature.value.level;
        points = shapePoints(x, y, width, height, feature.value.strand, feature.value.type);
        let invisible = false;
        if (feature.value.type === 'gene_feature') {
            rotateAngle = feature.value.rotate;
            if (isShowGeneLabel) {
                // Element width is too small and hide label at the edges
                if (width < 10 || start[0] >= rightCoord || end[0] <= leftCoord) {
                    invisible = true;
                } else {
                    invisible = false;
                }
            } else {
                invisible = true;
            }
            shape = graphic.clipPointsByRect(points, {
                x: params.coordSys.x,
                y: params.coordSys.y,
                width: params.coordSys.width,
                height: params.coordSys.height,
            });
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
                        text: feature.name,
                        fill: feature.itemStyle.color,
                        fontStyle: "normal",
                        fontSize: 12,
                        fontWeight: "bolder",
                    },
                },
                textConfig: {
                    position: "top",
                    distance: 18,
                    rotation: rotateAngle,
                    offset: "center",
                    local: true,
                },
            };
        } else if (feature.value.type === 'amplicon_feature') {
            shape = graphic.clipPointsByRect(points, {
                x: params.coordSys.x,
                y: params.coordSys.y,
                width: params.coordSys.width,
                height: params.coordSys.height,
            });
            return {
                type: "polygon",
                shape: {
                    points: shape,
                },
                style: api.style(),
                textContent: {},
                textConfig: {},
                invisible: !amplicon
            };
        } else if (feature.value.type === 'segment_feature') {
            if (width < 10 || start[0] >= rightCoord || end[0] <= leftCoord) {
                invisible = true;
            } else {
                invisible = false;
            }
            shape = graphic.clipPointsByRect(points, {
                x: params.coordSys.x,
                y: params.coordSys.y,
                width: params.coordSys.width,
                height: params.coordSys.height,
            });
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
                        text: feature.name,
                        fill: feature.itemStyle.color,
                        fontStyle: "normal",
                        fontSize: 12,
                        fontWeight: "bolder",
                    },
                },
                textConfig: {
                    position: "top",
                    distance: 18,
                    rotation: feature.value.rotate,
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

export {getGeneFeatureRenderer};