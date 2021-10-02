function selectDefaultSamples(samples){
    // Set default samples display 
    $("#selectedsamples").select2();
    $("#selectedsamples").val(samples);
    $("#scale").select2();
    $("#scale").val('log');
    document.getElementById('ymax').value = 100000
}

$(document).ready(function() {
    $('#selectedsamples').select2();
    $("#selectedsamples").on("change", function (e) {
        var selectData = $("#selectedsamples").select2("data");
        var samples_list = []
        for (var [key, entries] of selectData.entries()){
           samples_list.push(selectData[key].text)
        }
        updateOption(samples_list)
    });
});

$("#selectedsamples").select2({
    tags: true
});

$("#selectedsamples").on("select2:#selectedsamples", function (evt) {
    var element = evt.params.data.element;
    var $element = $(element);
    
    $element.detach();
    $(this).append($element);
    $(this).trigger("change");
});
