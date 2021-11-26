import {defaultTo} from "lodash/util";
const join = require('lodash/join');
const map = require('lodash/map');

/**
 * Define properties for gene/amplicon feature plot which is in the last index of grid
 * @type {{max_grid_height: number, rec_items_height: number, grid_height: string}}
 */
export const geneFeaturePlotProperties = {
    'max_grid_height': 80,
    'rec_items_height': 12,
    'grid_height': "15%"
}

/**
 * Dict of color for coloring variant position in the chart
 * @type {Dict[]}
 */
export const ntColor = {
    "A": "#ea5e48",
    "C": "#eaca48",
    "G": "#6ad82b",
    "T": "#2b87d8",
}

/**
 * Write tooltip information to HTML table
 * @param {string} headers - Header of table
 * @param {Array<Array<string>>} rows - Rows of table
 * @param {string} classes - Classes defined for table
 * @returns {string}
 */
function toTableHtml(headers, rows, classes) {
    var classTable = defaultTo(classes, "table");
    var out = '<table class="' + classTable + '"><thead>';
    out += join(
        map(headers, function (x) {
            return "<strong>" + x + "</strong>";
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
                        return "<td " + (i === 0 ? 'scope="row"' : "") + ">" + x + "</td>";
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

export {toTableHtml};