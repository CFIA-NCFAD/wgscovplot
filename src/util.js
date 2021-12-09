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
 * Function get Variant Comparison across samples
 * @param {Array<string>} samples - An array of samples name
 * @param {Array<Array<number>>} depths - Array of depths
 * @param {Array<Array<Object>>} variants - The dict of variants data
 * @param {number} position - Variant position
 * @param {string} currentSample - selected sample
 * @returns <Array<Array<string>> - Variant comparison across samples
 */
function getVariantComparison(samples, variants, depths, position, currentSample=""){
    var rows = [];
    var variantArr = [];
    for (var [i, element] of variants.entries()) {
        if (element.length){
            var isPOSExist = false;
            Object.values(element). forEach(values => {
                if (values['POS'] === position){
                    isPOSExist = true;
                    variantArr.push(values);
                }
            })
            if (!isPOSExist){
                variantArr.push({"sample": samples[i], "POS":position}); // sample has variant infor but no variant infor at this position
            }
        }
        else{
            variantArr.push({"sample": samples[i], "POS":position}); // sample has no variant information
        }
    }
    var unionKeys = [...new Set(variantArr.reduce((r, e) => [...r, ...Object.keys(e)], []))];
    unionKeys.push("Coverage Depth") // Add Coverage Depth row
    unionKeys.forEach(key => {
        var row = [];
        row.push(key);
        if (key === "Coverage Depth"){
            for (var [i, depthEle] of depths.entries()){
                row.push(depths[i][position-1].toLocaleString())
            }
        }else{
            variantArr.forEach(element => {
                if (element[key] !== undefined && element[key] !== null){
                    if (key === 'sample' && element[key] === currentSample) // Bold highlight selected sample
                        row.push(element[key].bold());
                    else{
                        row.push(element[key]);
                    }
                }else {
                    row.push("");
                }
            })
        }
        rows.push(...[row]);
    })
    return rows;
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

export {toTableHtml, getVariantComparison};