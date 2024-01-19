import {ECFeature, ECFormatterFeature, WgsCovPlotDB} from "../db";
import {getTextWidth, hexToHSL, shapePoints} from "../util";
import {graphic} from "echarts/core";
import {get, isEmpty, isNil} from "lodash";

interface EChartsParams {
  dataIndex: number;
  coordSys: {
    x: number;
    y: number;
    width: number;
    height: number;
  };
}

interface EChartsRenderItemAPI {
  // eslint-disable-next-line
  coord: (arg0: any[]) => [any, any] | [any];
  // eslint-disable-next-line
  style: () => any;
}

function getGeneFeatureRenderer(db: WgsCovPlotDB) {
  let yStart = 0;

  function renderGeneFeatures(params: EChartsParams, api: EChartsRenderItemAPI) {
    const {
      dataIndex,
      coordSys,
    } = params;
    const feature: ECFeature | undefined = get(db, ["echart_features", dataIndex]);
    if (isNil(feature)) {
      return null;
    }
    const leftCoord = coordSys.x;
    const rightCoord = coordSys.width + coordSys.x;
    const [startX, startY] = api.coord([feature.value.start, dataIndex]);
    if (dataIndex === 0) {
      yStart = startY;
    }
    const [endX] = api.coord([feature.value.end, dataIndex]);
    const height = db.chartOptions.geneLabelTextSize + 6;
    const width = endX - startX;
    const y = yStart - height / 2 - feature.value.level;
    const points = shapePoints(startX, y, width, height, feature.value.strand, feature.value.type);
    if (points === null) {
      return null;
    }
    const textWidth = getTextWidth(feature.name, `normal ${db.chartOptions.geneLabelTextSize}px Arial`);
    const shape = graphic.clipPointsByRect(points, coordSys);
    let invisible = false;
    if (feature.value.type === "gene") {
      if (db.chartOptions.showGeneLabels) {
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
            text: feature.name,
            fill: width > textWidth ? (hexToHSL(feature.itemStyle.color).l > 50 ? "black" : "white") : "black",
            fontStyle: "normal",
            fontSize: db.chartOptions.geneLabelTextSize,
            fontWeight: "normal",
            stroke: feature.itemStyle.color,
            lineWidth: 1,
          },
        },
        textConfig: {
          position: width <= textWidth ? [0, -db.chartOptions.geneLabelTextSize / 2 - 1] : "inside",
          rotation: width <= textWidth ? 0.5 : 0,
          local: false,
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
        invisible: !db.show_amplicons
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
            text: feature.name,
            stroke: feature.itemStyle.color,
            lineWidth: 1,
            fontStyle: "normal",
            fontSize: 16,
            fontWeight: "bolder",
          },
        },
        textConfig: {
          position: "inside",
          distance: 0,
          rotation: db.chartOptions.geneLabelRotation,
          origin: "center",
          local: true,
        },
      };
    } else {
      return null;
    }
  }

  return renderGeneFeatures;
}

export const getGeneFeatureSeries = (db: WgsCovPlotDB) => {
  const index = db.chartOptions.selectedSamples.length;
  return {
    type: "custom",
    animation: false,
    silent: false,
    xAxisIndex: index,
    yAxisIndex: index,
    renderItem: getGeneFeatureRenderer(db),
    labelLayout: {
      hideOverlap: false,
    },
    data: db.echart_features,
    tooltip: {
      trigger: "item",
      show: true,
      enterable: true,
      appendToBody: true,
      triggerOn: "mousemove",
      renderMode: "html",
      borderRadius: 6,
      borderWidth: 2,
      showContent: "true",
      position: "top",
      textStyle: {
        fontSize: 14,
        fontWeight: "normal",
      },
      formatter: function (feature: ECFormatterFeature) {
        if (isNil(db.segments)) {
          return <div class="w-full h-full">
            <p class="text-sm font-bold">{/*@once*/ feature.name}</p>
            <p class="text-sm">
              {/*@once*/ feature.value.start.toLocaleString()} - {/*@once*/ feature.value.end.toLocaleString()}
            </p>
            <p class="text-sm">{/*@once*/ (feature.value.end - feature.value.start + 1).toLocaleString()} bp</p>
            <button class="hover:ring py-1 px-2 rounded mt-2"
                    style={{
                      "background-color": /*@once*/ feature.color !== undefined ? feature.color : "#2a2a2a",
                      "color": /*@once*/ feature.color !== undefined ? (hexToHSL(feature.color).l > 50 ? "black" : "white") : "#fff",
                    }}
                    onClick={(e) => {
                      if (isNil(db.ref_seq)) {
                        return;
                      }
                      const seq = db.ref_seq.slice(feature.value.start - 1, feature.value.end);
                      // TODO: add reference ID/name to header; need to make sure that it is exported from the backend
                      const header = `>REFID|${feature.value.start}-${feature.value.end} ${feature.name}`;
                      navigator.clipboard.writeText(`${header}\n${seq}\n`);
                      // Indicate that the text has been copied
                      // Using a Solid signal will produce the following warning:
                      // "computations created outside a `createRoot` or `render` will never be disposed"
                      // so the button text is updated directly
                      e.currentTarget.innerHTML = "Copied seq!";
                    }}>
              Copy seq
            </button>
          </div>;
        }
      },
    },
  };
}

export function getRegionAmpliconDepthRenderer(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  function renderRegionAmpliconDepth({coordSys}: {coordSys: any}, api: any) {
    const [startX, startY] = api.coord([api.value(0), api.value(2)]);
    const [endX, endY] = api.coord([api.value(1), 1]);
    const rectShape = graphic.clipRectByRect(
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
      invisible: !db.show_amplicons
    };
  }

  return renderRegionAmpliconDepth;
}

export function getRegionAmpliconDepthSeries(db: WgsCovPlotDB) {
  // eslint-disable-next-line
  const ampliconDepthSeries: any[] = [];
  if (isEmpty(db.amplicon_depths))
    return ampliconDepthSeries;
  for (const [i, sample] of db.chartOptions.selectedSamples.entries()) {
    ampliconDepthSeries.push({
      type: "custom",
      xAxisIndex: i,
      yAxisIndex: i,
      renderItem: getRegionAmpliconDepthRenderer(db),
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
      data: db.amplicon_depths[sample],
    });
  }
  return ampliconDepthSeries;
}