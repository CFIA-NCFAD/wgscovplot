import {shapePoints} from "./shapePoints";
import {graphic} from "echarts/core";

/**
 * A closure is used for initial y_start for rendering gene/amplicon features
 * @type {function(): number}
 */
var yStart = (function (){
    var start;
    return function (){
        return start;
    };
})();

/**
 * Renders the gene/amplicon features shape
 * @param {boolean} isShowGeneLabel - Whether to show gene label or not
 * @param {Array<Dict[]>} geneFeatureAmpliconData - Array of dictionary geneFeature or amplicon data
 * @returns The shape of gene feature and/or amplicon feature
 */
function getGeneFeatureRenderer(isShowGeneLabel, geneFeatureAmpliconData) {
    /**
     *
     * @param params - Echarts params
     * @param api - Echarts api
     */
    function renderGeneFeatures(params, api) {
        var points, shape, rotateAngle;
        var start, end, height, width, x, y;
        var categoryIndex = params.dataIndex;
        var feature = geneFeatureAmpliconData[categoryIndex];
        start = api.coord([feature.value.start, categoryIndex]);
        if (categoryIndex === 0) {
            yStart = start[1];
        }
        end = api.coord([feature.value.end, categoryIndex]);
        height = geneFeatureProperties["rec_items_height"];
        width = end[0] - start[0];
        x = start[0];
        y = yStart - height / 2 - feature.value.level;
        points = shapePoints(x, y, width, height, feature.value.strand, feature.value.type);
        if (feature.value.type === 'gene_feature') {
            rotateAngle = (feature.value.strand === 1) ? 0.7 : -0.7;
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
                style: api.style({}),
                textContent: {
                    type: "text",
                    invisible: !isShowGeneLabel,
                    style: {
                        text: feature.name,
                        fill: feature.itemStyle.color,
                        fontStyle: "normal",
                        fontSize: 10,
                        fontWeight: "bolder",
                    },
                },
                textConfig: {
                    position: "top",
                    distance: 20,
                    rotation: rotateAngle,
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
            };
        } else {
            return null;
        }
    }
    return renderGeneFeatures;
}

export {getGeneFeatureRenderer};