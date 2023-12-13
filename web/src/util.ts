import {join, map} from "lodash";

/**
 * Define properties for gene/amplicon feature plot which is in the last index of grid
 */
export const FEATURE_PLOT_PROPS = {
  "max_grid_height": 80,
  "rec_items_height": 12,
  "grid_height": "15%"
};

/**
 * Define color for segmented gene segments
 */
export const FLU_SEGMENT_COLOURS = {
  "1_PB2": "#A6CEE3",
  "2_PB1": "#1F78B4",
  "3_PA": "#B2DF8A",
  "4_HA": "#33A02C",
  "5_NP": "#FB9A99",
  "6_NA": "#E31A1C",
  "7_M": "#FDBF6F",
  "8_NS": "#FF7F00"
};

/**
 * Write tooltip information to HTML table
 * @param {string[]} headers - Header of table
 * @param {Array<Array<string>>} rows - Rows of table
 * @param {string} classes - Classes defined for table
 * @returns {string}
 */
export const toTableHtml = ({headers, rows, classes = "table small"}
                              : { headers: string[], rows: string[][], classes?: string }) => {
  let out = `<table class="${classes}"><thead>`;
  out += join(
    map(headers, function (x) {
      return `<strong>${x}</strong>`;
    }),
    ""
  );
  out += "</thead><tbody>";
  out += join(
    map(rows, function (xs) {
      return (
        "<tr>" +
        join(
          map(xs, function (x, i) {
            return `<td ${i === 0 ? "scope=\"row\"" : ""}><samp>${x}</samp></td>`;
          }),
          ""
        ) +
        "</tr>"
      );
    }),
    ""
  );
  out += "</tbody></table>";
  return out;
}

/**
 * Get the regions in which depth < threshold
 */
export function getCoordsInterval(depths: number[], threshold: number) {
  let foundInterval = false;
  let firstCoord = 0;
  let n = depths.length;
  const out = [];
  let i = 0;
  for (; i < n; ++i) {
    if (depths[i] < threshold) {
      if (!foundInterval) {
        firstCoord = i;
      }
      foundInterval = true;
    } else {
      if (foundInterval) {
        out.push([firstCoord + 1, i]);
      }
      foundInterval = false;
    }
  }
  // add the last interval
  if (foundInterval) {
    out.push([firstCoord + 1, i]);
  }
  return out;
}


/**
 * Convert hex color to HSL
 *
 * Example:
 *
 *  hexToHSL('#FF0000') // {h: 0, s: 1, l: 0.5}
 *
 * courtesy of https://stackoverflow.com/a/62390801
 */
export function hexToHSL(hex: string): { h: number; s: number; l: number; } {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (result === null || result.length !== 4) {
    return {h: 0, s: 0, l: 0};
  }
  let r = parseInt(result[1], 16);
  let g = parseInt(result[2], 16);
  let b = parseInt(result[3], 16);

  r /= 255
  g /= 255
  b /= 255;
  const max = Math.max(r, g, b)
  const min = Math.min(r, g, b);
  let h, s, l = (max + min) / 2;

  if (max == min) {
    h = s = 0; // achromatic
  } else {
    const d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;
      case g:
        h = (b - r) / d + 2;
        break;
      case b:
        h = (r - g) / d + 4;
        break;
    }
    if (h !== undefined) {
      h /= 6;
    }
  }

  s = s * 100;
  s = Math.round(s);
  l = l * 100;
  l = Math.round(l);
  if (h !== undefined) {
    h = Math.round(360 * h);
  } else {
    h = 0;
  }
  return {h, s, l};
}

export function shapePoints(
  x: number,
  y: number,
  width: number,
  height: number,
  strand: number,
  type: string
): number[][] | null {
  let taperOffset: number;
  // Element width is too small
  if (width < 10) {
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


export function getTextWidth(text: string, font: string): number {
  // re-use canvas object for better performance
  const canvas = document.createElement("canvas");
  const context = canvas.getContext("2d");
  if (context === null) {
    return 0;
  }
  context.font = font;
  const metrics = context.measureText(text);
  return metrics.width;
}
