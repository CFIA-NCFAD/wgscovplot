/**
 * Write tooltip information to HTML table
 * @param {string} headers - Header of table
 * @param {Array<Array<string>>} rows - Rows of table
 * @param {string} classes - Classes defined for table
 * @returns {string}
 */
function toTableHtml(headers, rows, classes) {
    var classTable = _.defaultTo(classes, "table");
    var out = '<table class="' + classTable + '"><thead>';
    out += _.join(
        _.map(headers, function (x) {
            return "<strong>" + x + "</strong>";
        }),
        ""
    );
    out += "</thead><tbody>";
    out += _.join(
        _.map(rows, function (xs) {
            return (
                "<tr>" +
                _.join(
                    _.map(xs, function (x, i) {
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