import {defaultTo} from "lodash/util";
const join = require('lodash/join');
const map = require('lodash/map');

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