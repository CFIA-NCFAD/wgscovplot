function renderenv(){
    var render_env = document.getElementById('renderenv').value;
    var isChecked=document.getElementById("toggledarkmode").checked;
    var samples = [];
    
    for (const [key, entries] of Object.entries(window.samples)){
        if (key < default_num_chart){
            samples.push(entries)
        }
    }
    var mode =  'white';
    if (isChecked)
        mode = 'dark'
    else
        mode = 'white'
    if (render_env==='canvas'){
        echarts.dispose(chart) // destroy chart instance and re-init chart
        $chart = document.getElementById('chart');
        chart = echarts.init($chart, mode , {renderer: 'canvas'});
        chart.setOption(option = getOption())
        selectDefaultSamples(samples)
        chartDispatchAction()
    }
    else{
        echarts.dispose(chart)
        $chart = document.getElementById('chart');
        chart = echarts.init($chart, mode , {renderer: 'svg'});
        chart.setOption(option = getOption())
        selectDefaultSamples(samples)
        chartDispatchAction()
    }
}

$('#toggledarkmode').change(function() {
    var render_env = 'canvas'
    if (document.getElementById('renderenv'))
        render_env = document.getElementById('renderenv').value;
    var samples = [];
    for (const [key, entries] of Object.entries(window.samples)){
        if (key < default_num_chart){
            samples.push(entries)
        }
    }
    if ($(this).prop('checked')){
        echarts.dispose(chart) // destroy chart instance and re-init chart
        chart = echarts.init($chart, 'dark', {renderer: render_env});
        chart.setOption(option = getOption())
        selectDefaultSamples(samples)
        chartDispatchAction()
    }
    else{

        echarts.dispose(chart) // destroy chart instance and re-init chart
        chart = echarts.init($chart, 'white', {renderer: render_env});
        chart.setOption(option = getOption())
        selectDefaultSamples(samples)
        chartDispatchAction()
    } 
})


function chartDispatchAction(){
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
}


function updateChartHeight(val) {
    document.getElementById('chartheight').value=val +'%'; 
    var grid_option = chart.getOption().grid;
    for (var i = 0; i < grid_option.length; i++){
        if (i < grid_option.length - 1) // don't change height of DNA feature charts
            grid_option[i]['height'] = val +'%'
        if (i > 0){
            grid_option[i]['top'] = parseInt(grid_option[i -1 ]['top'].replace('%','')) + parseInt(grid_option[i - 1]['height'].replace('%','')) + 4 + '%'
        }
    }
    chart.setOption({grid:grid_option})
    
}

function updateChartLeft(val) {
    document.getElementById('chartleft').value=val+'%'; 
    var grid_option = chart.getOption().grid;
    for (var i = 0; i < grid_option.length; i++){
        grid_option[i]['left'] = val +'%'
    }
    chart.setOption({grid:grid_option})
    
}


function updateChartRight(val) {
    document.getElementById('chartright').value=val+'%'; 
    var grid_option = chart.getOption().grid;
    for (var i = 0; i < grid_option.length; i++){
        grid_option[i]['right'] = val +'%'
    }
    chart.setOption({grid:grid_option})
    
}

function updateChartTop(val) {
    document.getElementById('charttop').value=val+'%'; 
    var grid_option = chart.getOption().grid;
    for (var i = 1; i < grid_option.length; i++){
        grid_option[i]['top'] = parseInt(grid_option[i-1]['top'].replace('%','')) + parseInt(grid_option[i-1]['height'].replace('%','')) + parseInt(val) + '%'
    }
    chart.setOption({grid:grid_option})
    
}

