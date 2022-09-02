/**
 * Get toolbox option
 * @returns {Object}
 */
function getToolbox() {
    let toolBox;
    toolBox = {
        show: "true",
        feature: {
            dataView: {
                readOnly: false,
            },
            restore: {},
            saveAsImage: {
                name: "wgscovplot",
            },
        },
    };
    return toolBox;
}

export {getToolbox};