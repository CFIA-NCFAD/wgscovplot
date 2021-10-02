function setScale(){
    var scaletype = document.getElementById('scale').value;
    var yaxis_option = chart.getOption().yAxis;
    if (scaletype == 'value'){
        for(let i = 0 ; i < grid_length; i++){
            yaxis_option[i].type = scaletype
            yaxis_option[i].min = 0
            yaxis_option[i].max = 40000
        }
        chart.setOption({yAxis: yaxis_option})
        document.getElementById('ymax').value = 40000
    }
    else{
        for(let i = 0 ; i < grid_length; i++){
            yaxis_option[i].type = scaletype
            yaxis_option[i].min = 1
            yaxis_option[i].max = 100000
        }
        chart.setOption({yAxis: yaxis_option})
        document.getElementById('ymax').value = 100000
    }
}

function setYMax(){
    var ymax = document.getElementById('ymax').value;
    var yaxis_option = chart.getOption().yAxis;
    for(let i = 0 ; i < grid_length; i++){
        yaxis_option[i].max = ymax
    }
    chart.setOption({yAxis: yaxis_option})
}

