$('#togglegenelabel').change(function() {
    var series_option = chart.getOption().series;
    gene_features_grid_index = series_option.length - 1
    if ($(this).prop('checked')){
        series_option[gene_features_grid_index]['renderItem'] = renderGeneFeatures
        chart.setOption({series: [...series_option]})
    }
    else{
        series_option[gene_features_grid_index]['renderItem'] = updateGeneFeatures
        chart.setOption({series: [...series_option]})
    } 
})

function updateGeneFeatureHeight(val) {
    document.getElementById('genefeatureheight').value=val +'%'; 
    var grid_option = chart.getOption().grid;
    gene_features_grid_index = grid_option.length - 1
    grid_option[gene_features_grid_index]['height'] = val +'%'
    chart.setOption({grid:grid_option})
}
