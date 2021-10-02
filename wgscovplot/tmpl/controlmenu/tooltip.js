$('#toggletooltip').change(function() {
    if ($(this).prop('checked')){
        chart.setOption({
            tooltip:{showContent: true}
        });
    }
    else{
        chart.setOption({
            tooltip:{showContent: false}
        });
    }
})
