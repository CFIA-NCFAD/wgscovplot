$('#toggleslider').change(function() {
    if ($(this).prop('checked')){
        chart.setOption({
            dataZoom: [
            {
                type: "inside",
                filterMode: 'none',
                xAxisIndex: [...Array(grid_length + 1).keys()]
            },
            {
                show: true,
                filterMode: 'none',
                xAxisIndex: [...Array(grid_length + 1).keys()],
                type: "slider"
            }
        ]
        });
    }
    else{
        chart.setOption({
            dataZoom: [
            {
                type: "inside",
                filterMode: 'none',
                xAxisIndex: [...Array(grid_length + 1).keys()]
            },
            {
                show: false,
                type: "slider"
            }
        ]
        });
    } 
})

function setDataZoom(){
    var start = document.getElementById('start_pos').value
    var end = document.getElementById('end_pos').value
    chart.dispatchAction({
        type:'dataZoom',
        startValue: start,
        endValue: end
    })
}

function resetDataZoom(){
    document.getElementById('start_pos').value = 1 
    document.getElementById('end_pos').value = ref_len
    chart.dispatchAction({
        type:'dataZoom',
        start: 0,
        end: 100
    })
}


chart.on('click', function(params){
    document.getElementById('start_pos').value = params.value[1]
    document.getElementById('end_pos').value = params.value[2]
    chart.dispatchAction({
        type:'dataZoom',
        startValue: params.value[1],
        endValue: params.value[2]
    })
})

chart.on('dblclick', function(params){
    document.getElementById('start_pos').value = 1
    document.getElementById('end_pos').value = ref_len
    chart.dispatchAction({
        type:'dataZoom',
        start: 0,
        end: 100
    })
})
