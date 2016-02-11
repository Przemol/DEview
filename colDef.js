[
{ "title": 'Gene ID',
  "render": function ( data, type, row ) {
      return '<a href="' + $('#annosoeurce').val() + data + '" target="_new">' + data + '</a>'
  }
},
{ "title": 'Gene' },
{ "title": 'Seq ID' },
{ "title": 'Type' },
  { "title": '',
    "className": "table-center",
    "width": "1%",
    "orderable": false,
    "searchable": false,
    "render": function ( data, type, row ) {
      if (data.length) {
        return '<button class="btn btn-info btn-xs no_select" rel="popover" data-content="'+ data +'" data-original-title="'+row[1] + ' (' +row[0]+')"'+ 
          'onmouseover="$(this).popover({trigger: \'hover\', placement: \'left\'}).popover(\'show\')">'+'<span class="glyphicon glyphicon glyphicon-info-sign" aria-hidden="true"></span></button>'
      } else {
        return '<div style="width: 29px"></div>'
      }
    }
  },


{ "title": 'Base mean', type: "number-range" },
{ "title": 'MAP FC' },
{ "title": 'LogFC' },
{ "title": 'FC error' },
{ "title": 'Statistic val' },
{ "title": 'p-value',
  "render": function ( data, type, row ) {
    return '<span style="white-space: nowrap;">' + Number(data).toPrecision(2) + '</span>';
  }
},
{ "title": 'FDR',
  "render": function ( data, type, row ) {
    return '<span style="white-space: nowrap;">' + Number(data).toPrecision(2) + '</span>';
  }
},

{ "title": 'Plot',
  "className": "table-center",
  "orderable": false,
  "searchable": false,
  "width": "1%",
  "render": function ( data, type, row ) {
      return '<button class="btn btn-xs btn-success no_select"'+
      'onClick="$(\'[data-value=Plot]\').click(); Shiny.shinyapp.sendInput({\'plot\':\''+row[0]+'\'});'+
      '"><span class="glyphicon glyphicon-stats" aria-hidden="true"></span></button>'
  }
}

]
