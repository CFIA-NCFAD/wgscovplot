import {getFluGrids} from "./getFluGrids";

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

function getYAxisName(samples, segments, samplesIndex, segmentsIndex){

   let yAxisName = '';
   let horizontalOffset = 4.0;
   let padLeft = 4.0;
   let samplePlotWidth = (100.0 - padLeft)/segments.length - horizontalOffset;
   if (samplePlotWidth > 10){
       yAxisName = samples[samplesIndex] + ' - Segment ' + segments[segmentsIndex];
   }else{
       yAxisName =  samples[samplesIndex].substring(0,10) + '...';
   }
   return yAxisName;
}

function getFluYAxes(samples, segments) {
    let axes = [];
    let gridIndex = 0;
    for (let i = 0; i < samples.length; i++) {
        for (let j = 0; j < segments.length; j++) {
            axes.push({
                type: 'log',
                gridIndex: gridIndex,
                name: getYAxisName(samples, segments, i, j),
                nameTextStyle: {
                    fontStyle: "normal",
                    fontWeight: "normal",
                    fontSize: 10,
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