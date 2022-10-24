/**
 * Get toolbox option
 * @returns {Object}
 */
function getToolbox() {
    let toolBox;
    toolBox = {
        show: "true",
        feature: {
            saveAsImage: {
                name: "wgscovplot",
            },
            restore: {},
            dataView: {
                readOnly: false,
            }
        },
    };
    return toolBox;
}

export {getToolbox};