$(function() {
    $('#downloadDataFlt').click(function(){ 
        var zz = $($('#data').find('table')[1]).DataTable().ajax.params();
        zz.length = -1; 
        Shiny.shinyapp.sendInput({sel:zz});
    })
});