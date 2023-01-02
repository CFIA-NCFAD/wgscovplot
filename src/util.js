import {map} from "lodash/collection";
import {join} from "lodash/array";
let isBase64 = require("is-base64");

/**
 * Define properties for gene/amplicon feature plot which is in the last index of grid
 */
export const FEATURE_PLOT_PROPS = {
    "max_grid_height": 80,
    "rec_items_height": 12,
    "grid_height": "15%"
};

/**
 * Dict of color for coloring variant position in the chart
 */
export const NT_COLOURS = {
    "A": "#ea5e48",
    "C": "#eaca48",
    "G": "#6ad82b",
    "T": "#2b87d8",
};

/**
 * Define color for flu gene segments
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

//Now support both base64 and uncompress array
function toFloat32Array(b64String) {
    if (isBase64(b64String)) {
        let f32a = new Float32Array(new Uint8Array([...window.atob(b64String)].map(c => c.charCodeAt(0))).buffer);
        let depthArr = Array.from(f32a);
        return depthArr;
    }
    return b64String;
}


/**
 * Write tooltip information to HTML table
 * @param {string[]} headers - Header of table
 * @param {Array<Array<string>>} rows - Rows of table
 * @param {string} classes - Classes defined for table
 * @returns {string}
 */
function toTableHtml(
    {
        headers,
        rows,
        classes = "table small",
    }) {
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
 *
 * @param {Array<number>} depths
 * @param {number} threshold
 */
function getCoordsInterval(depths, threshold) {
    let coords = [];
    let foundInterval = false;
    let firstCoord = 0;
    let count = 0;
    for (let i = 0; i < depths.length; i++) {
        if (depths[i] < threshold) {
            firstCoord = i - count;
            count += 1;
            foundInterval = true;
            if (i === depths.length - 1) {
                let end = firstCoord + count;
                let start = firstCoord + 1;
                coords.push({start, end}); // pos in index 1
            }
        } else {
            if (foundInterval) {
                let start = firstCoord + 1;
                let end = firstCoord + count;
                coords.push({start, end}); // pos in index 1
            }
            foundInterval = false;
            count = 0;
        }
    }
    return coords;
}

export {toTableHtml, getCoordsInterval, toFloat32Array};