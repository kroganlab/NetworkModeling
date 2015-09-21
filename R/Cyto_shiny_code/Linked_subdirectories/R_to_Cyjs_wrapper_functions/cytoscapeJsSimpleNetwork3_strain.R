cytoscapeJsSimpleNetwork2 <- function(nodeEntries, edgeEntries, 
         standAlone=FALSE, layout="cose", data_shape = 'ellipse', 
         #                                    height=600, width=600, injectCode="") {
#         height=700, width=800, injectCode="") {
  height=900, width=1100, injectCode="") {
    # Create webpage
    PageHeader <- "
    <!DOCTYPE html>
    <html>
    <head>
    <meta name='description' content='[An example of getting started with Cytoscape.js]' />
    
 //   <script src='http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js'></script>
 //   <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js'></script>
 //   <script src='file:cytoscape.min.js'></script>
 //    <script src='http://localhost/cytoscape.min.js'></script>
      <script src='cytoscape.min.js'></script>
 //   <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/arbor.js'></script>
//    <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cola.v3.min.js'></script>
//    <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/springy.js'></script>
//    <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/dagre.js'></script>
 //       <script src='http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/arbor.js'></script>
    
    <meta charset='utf-8' />
    <title> Influenza </title>"
    
    if(standAlone) {
        NetworkCSS <- "<style>
        #cy {
        height: 100%; 
        width: 100%;
        position: absolute; 
        left: 0;
        top: 200;
        border: 2px solid; 
    }
        </style>"    
} else {
    NetworkCSS <- paste0("<style>
                         #cy {
                         height: ", height, "px; 
                         width: ", width, "px;
                         position: relative;
                         left: 0;
                         top: 200;
                         border: 2px solid; 
}</style>")        
    }

# Main script for creating the graph
MainScript <- paste0("
                     <script>
                     $(function(){ // on dom ready	
                     $('#cy').cytoscape({
                    boxSelectionEnabled: true,
                    liveUpdate : false,
                    hideLabelsOnViewport: true,
                     motionBlur: true,
                     style: cytoscape.stylesheet()
                     .selector('node')
                     .css({
 //       'width': '80px',
         'width': 'mapData(log2fc_attrib, 0, 10.0, 0, 100.0)',
 //       'height': '80px',
         'height': 'mapData(log2fc_attrib, 0, 10.0, 0, 100.0)',
        'content': 'data(name)',
        'pie-size': '85%',
 //        'pie-size': 'mapData(log2fc_attrib, 0, 10.0, 0, 100.0)',
        'background-color': '#FFFFFF',
        'background-color': 'data(bkgrnd_highlight)',
 //       'pie-1-background-color': '#E8747C',
      'pie-1-background-color': '#FF0000',
        'pie-1-background-opacity': 'mapData(h1n1_time, 0, 1.0, 0, 1.0)',
        'pie-1-background-size': 'mapData(H1N1, 0, 18, 0, 100)',
 //       'pie-2-background-color': '#74CBE8',
 //     'pie-2-background-color': '#0000FF',
        'pie-2-background-color': '#000099',
        'pie-2-background-opacity': 'mapData(h3n2_time, 0, 1.0, 0, 1.0)',
        'pie-2-background-size': 'mapData(H3N2, 0, 18, 0, 100)',
 //       'pie-3-background-color': '#74E883',
       'pie-3-background-color': '#006600',
        'pie-3-background-opacity': 'mapData(h5n1_time, 0, 1.0, 0, 1.0)',
        'pie-3-background-size': 'mapData(H5N1, 0, 18, 0, 100)',
   //             'font-size': 16,
                'font-size': 'mapData(log2fc_attrib, 0, 10.0, 0, 24.0)',
                    'text-valign': 'center',
                     'color': 'black',
                     'text-outline-width': 1,
  //                    'shape': 'hexagon',
                    'shape':  '", data_shape, "',
                      'border-color': 'black',
                      'border-width': 1,
                      'border-style': 'solid',
 //                    'text-outline-color': 'white'
                   'text-outline-color': 'data(sign_color)',
                     })
                     .selector('edge')
                     .css({
                     'line-color': 'data(color)',
//                    'line-color': 'black',
                        'width': 1,
                        'line - width': 1,
                     'curve-style': 'haystack',
//                     'source-arrow-color': 'data(color)',
//                     'target-arrow-color': 'data(color)',
//                     'source-arrow-shape': 'data(edgeSourceShape)',
 //                    'target-arrow-shape': 'data(edgeTargetShape)'
                     })
                     .selector(':selected')
                     .css({
                     'background-color': 'black',
                     'line-color': 'black',
 //                    'target-arrow-color': 'black',
//                     'source-arrow-color': 'black'
                     })
                     .selector('.faded')
                     .css({
                     'opacity': 0.25,
                     'text-opacity': 0
                     }),
                     
                     elements: {
                     nodes: [",
                     nodeEntries,
                     "],
                     edges: [",
                     edgeEntries,
                     "]
                     },
                     
                     layout: {
                     name: '", layout, "',
                     padding: 10
                     },
                     
                     ready: function() {
                     window.cy = this; 
                     
                     //Injected options
                     ",
                     injectCode 
                     , "
                     
                     cy.on('tap', 'node', function(){
                     if(this.data('href').length > 0) {
                     window.open(this.data('href'));
                     }
                     
                     //console.log(this.data('href'));
                     });
                     }
                     }); 
                     }); // on dom ready
                     </script>")

PageBody <- "</head><body><div id='cy'></div>"
PageFooter <- "</body></html>"	

if(standAlone) {
    results <- paste0(PageHeader, NetworkCSS, MainScript, PageBody, PageFooter)	        
    return(results)
} else {
    results <- paste0(NetworkCSS, MainScript)    
    cat(results)
}
}