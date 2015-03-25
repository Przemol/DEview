$(function() {
    $('#downloadDataFlt').click(function(){ 
        var zz = $('table').DataTable().ajax.params(); 
        zz.length = -1; 
        Shiny.shinyapp.sendInput({sel:zz});
    })
});