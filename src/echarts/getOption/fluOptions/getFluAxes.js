
function getFluXAxes(samples, segments, depths) {
    let axes = [];
    for (let i = 0; i < samples.length * segments.length; i++) {
        axes.push({
            type: "value",
            gridIndex: i,
            min: 1,
            max: depths[i].length,
            axisLabel: {
                show:true,
                interval: "auto",
                hideOverlap: true
            }
        });
    }
    return axes;
}
function getFluYAxes(samples, segments) {
    let axes = [];
    let gridIndex = 0;
    for (let i = 0; i < samples.length; i++) {
        for (let j = 0; j < segments.length; j++) {
            axes.push({
                type: 'log',
                gridIndex: gridIndex,
                name: samples[i] + ' - Segment ' + segments[j],
                nameTextStyle: {
                    fontStyle: "normal",
                    fontWeight: "normal",
                    fontSize: 13,
                    width: 50,
                    overflow: 'truncate',
                    ellipsis: '...'
                },
                nameLocation: "end",
                nameRotate: 0.01,
                min: 1,
                max: 100000,
                minorSplitLine: {
                    show: true,
                },
            });
            gridIndex = gridIndex + 1;
        }
    }
    return axes;
}

export {getFluXAxes, getFluYAxes};