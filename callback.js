$('table tfoot th').slice(0,2).each( function () {
    var title = $('table thead th').eq( $(this).index() ).text();
    var width = $('table thead th').eq( $(this).index() ).width()+25;
    $(this).html( '<input type=\"text\" placeholder=\"'+title+'\" style=\"width:'+width+'px;\" />' );
} );

$('table tfoot th').slice(3,9).each( function () {
    var title = $('table thead th').eq( $(this).index() ).text();
    var width = 50;
    $(this).html( '<input class=\"min\" type=\"text\" placeholder=\"'+'min'+'\" style=\"width:'+width+'px;\" /><br />' +
                  '<input class=\"max\" type=\"text\" placeholder=\"'+'max'+'\" style=\"width:'+width+'px;\" />' );
} );
$('table tfoot th').css('text-align', 'right');
$('table tfoot th').css('padding', 5);

table.columns().eq( 0 ).each( function ( colIdx ) {
    $( 'input', table.column( colIdx ).footer() ).on( 'keyup change', function () {
        if(this.className == 'min') {
            var flt = $(this).val() + ',' + $(this).siblings('.max').val();
            table.column( colIdx ).search( flt ).draw();
        } else if(this.className == 'max') {
            var flt = $(this).siblings('.min').val() + ',' + $(this).val();
            table.column( colIdx ).search( flt ).draw();
        } else {
            table.column( colIdx ).search( this.value ).draw();
        }
    } );
} );
            