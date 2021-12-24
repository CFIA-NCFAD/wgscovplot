import {toTableHtml} from "../../util";

function getTooltipHeatmap(samples, mutations, variants){
    return [
        {
            enterable: true,
            appendToBody: true,
            renderMode: "html",
            showContent: true,
            confine: true,
            formatter: function (params){
                var output = '';
                var mutationName = mutations[params.value[0]];
                var sample = samples[params.value[1]];
                var rows = [];
                var isExist = false;
                if (variants[sample] !== undefined && variants[sample] !== null){
                    Object.values(variants[sample]).forEach(element => {
                        if (element["mutation"] === mutationName) {
                            for (const [key, value] of Object.entries(element)) {
                                rows.push(...[[key, value]])
                            }
                            isExist = true;
                        }
                    })
                    if (isExist){
                        output += toTableHtml(["", ""], rows, "table table-hover table-bordered table-responsive-md");
                    }
                    else{
                        output += ("Not detected mutation " +  mutationName.bold()  + " in sample " + sample.bold())
                    }
                }
                else{
                    output += "Sample " + sample.bold() + "has no variant inforamtion";
                }
                return output
            }
        }
    ]
}


function getVariantHeatmapOption(samples, mutations, variantMatrix, variants) {
    var chartOptions = {
        xAxis: {
            type: "category",
            data: mutations,
            splitArea: {
                show: true
            },
            position: "bottom",
            axisLabel:{
                rotate: -45,
            }
        },
        yAxis: {
            type: "category",
            data: samples,
            splitArea:{
                show: true
            },
        },
        dataZoom: [
            {
                type: "inside"
            },
            {
                type: "slider"
            }
        ],
        visualMap: {
            min: 0,
            max: 1,
            calculable: true,
            orient: 'vertical',
            left: 'right',
            top: '5%',
        },
        series: [
            {
                type: "heatmap",
                data: variantMatrix,
                label: {
                    show: false
                },
                emphasis: {
                    itemStyle: {
                        shadowBlur: 10,
                        shadowColor: 'rgba(0, 0, 0, 0.5)'
                    }
                }
            }
        ],
        tooltip: getTooltipHeatmap(samples, mutations, variants),
        toolbox: {
            show: "true",
            feature: {
                dataView: {
                    readOnly: false,
                },
                restore: {},
                saveAsImage: {
                    name: "wgscovplot_variant_matrix",
                },
            },
        },
        grid: {
            top: "5%",
            left: "15%",
            height: "50%",
        }
    };
    return chartOptions;
}

export {getVariantHeatmapOption};