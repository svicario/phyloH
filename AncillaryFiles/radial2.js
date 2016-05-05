var exe = require("exelixis");
var root_chain_svg;
//var Zip = require("adm-zip");
var createTree = exe.createTree;

updateTree = exe.updateTree;

d3.select("#nRadius").on("input", function() {
  updateCircleRadius(+this.value);
});

d3.select("#save").on("click", function(){
  var html = d3.select("svg")
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;

  //console.log(html);
  var imgsrc = 'data:image/svg+xml;base64,'+ btoa(html);
  var img = '<img src="'+imgsrc+'">'; 
  d3.select("#svgdataurl").html(img);


  var canvas = document.querySelector("canvas"),
      context = canvas.getContext("2d");

  var image = new Image;
  image.src = imgsrc;
  image.onload = function() {
      context.drawImage(image, 0, 0);

      var canvasdata = canvas.toDataURL("image/png");

      var pngimg = '<img src="'+canvasdata+'">'; 
      d3.select("#pngdataurl").html(pngimg);

      var a = document.createElement("a");
      a.download = "sample.png";
      a.href = canvasdata;
      a.click();
  };

});

updateCircleRadius(10);

initHoverTooltip();

testRootChain();

//testContextMenu();

/* Try to change the opts in your console */

opts = {

        el : document.getElementById("yourDiv"),
        tree : {
            //data : "(homo_sapiens:12,(mus_musculus:12,(danio_rerio:13,(pan_troglodytes:9,(taeniopygia_guttata:10,callithrix_jacchus:11):12):12):10);",
            data : "();",  
            width : 2000,
            scale : true,
            layoutInput : "radial",
        },
        label : {
            fontsize : 12,
            usePics : true, 
            pics : {
                pictureWidth : 30,
                pictureHeight : 40,
            },
        },
        nodes : {
            toggle : true, //allows onClickEvent
            select: true, //allows selections
            size : 5,
            fill : "orange",
            stroke : "yellow",
            selectedFill : "steelblue",
            selectedSize : 4,
        },
        histograms: {
            full_length: 70
        },
        upstream: []
};

opts.root_chain = d3.select("#root_chain")
                    .append("svg")
                    .attr("width", window.innerWidth)
                    .attr("height", 100);

tree = createTree(opts);

//updateTree(tree,opts);
json_to_newick = function (json) {
    function nested(nest)
    {
        var subtree = "";

        if(nest.hasOwnProperty('children')){
            var children = [];
            nest.children.forEach(function(child){
                var subsubtree = nested(child);
                children.push(subsubtree);
            });
            var substring = children.join();
            if(nest.hasOwnProperty('name')){
                subtree = "("+substring+")" + nest.name;
            }
            if(nest.hasOwnProperty('branch_length')){
                subtree = subtree + ":"+nest.branch_length;
            }
        }
        else 
        {
            var leaf = "";
            if(nest.hasOwnProperty('name')){
                leaf = nest.name;
            }
            if(nest.hasOwnProperty('branch_length')){
                leaf = leaf + ":"+nest.branch_length;
            }
            subtree = subtree + leaf;
        }
        return subtree;
    }
    return nested(json) +";";
};

computeInnerLeavesHistograms = function (tree) {
    function nested(nest)
    {
        if(nest.hasOwnProperty('children')){
            var children = [];
            nest.children.forEach(function(child){
                nested(child);
            });

            var count = nest.children.length;
            if (count>0)
            {
                //debugger;
                opts.table_hist.histograms[nest.name] = [];
                for (var l in opts.table_hist.labels)
                    opts.table_hist.histograms[nest.name].push(0);
                for (var cidx = 0; cidx<count; cidx++)
                {
                    var el = nest.children[cidx];
                    for (var lidx=0;lidx<opts.table_hist.histograms[nest.name].length;lidx++)
                    {
                        opts.table_hist.histograms[nest.name][lidx]+=
                            opts.table_hist.histograms[el.name][lidx];    
                    }
                }
            }
        }
    }
    nested(tree);
    //debugger;
    console.log(tree);
};

function normalize_hist(table_hist)
{
    var accumulators = [];
    for (var l in table_hist.labels)
        accumulators.push(0.0);
    for (var key in table_hist.histograms) {
      if (table_hist.histograms.hasOwnProperty(key)) {
        var idx = 0;
        var hs = table_hist.histograms[key];
        for (var l in table_hist.labels)
        {
            accumulators[idx]+= hs[idx];
            idx += 1;
        }
      }
    }
    for (var key in table_hist.histograms) {
      if (table_hist.histograms.hasOwnProperty(key)) {
        var idx = 0;
        for (var l in table_hist.labels)
        {
            table_hist.histograms[key][idx] = table_hist.histograms[key][idx]/accumulators[idx];
            idx += 1;
        }
      }
    }
    return table_hist;
}


parseTableHist = function (lines) {
    var table_hist = {"histograms":{}};
    for(var line = 0; line < lines.length; line++)
    {
        //console.log(lines[line]);
        var l = lines[line];
        var parts = l.split(",");
        var label = parts[0];
        var idx = 0;
        switch(label) {
            case "LABELS":
                table_hist.labels = parts.slice(1,parts.length);
                break;
            case "COLORS":
                table_hist.colors = parts.slice(1,parts.length);
                break;
            case ",,,":
                break;
            default:
                if (label.length>0)
                {
                    // console.log(l);
                    // console.log(parts);
                    table_hist.histograms[label] = [];
                    for (idx=1;idx<parts.length;idx++)
                        table_hist.histograms[label].push(parseFloat(parts[idx]));
                }
        }
    }
    // console.log(table_hist);
    table_hist = normalize_hist(table_hist);
    // console.log('normalized');
    // console.log(table_hist);
    return table_hist;
};

// Sembra il nome di un medicinale...
parseXItol = function (lines) {
    var xitol = {};
    for(var line = 0; line < lines.length; line++)
    {
        // console.log(lines[line]);
        var l = lines[line];
        var parts = l.split("\t");
        var node_name = parts[0];
        var var_type = parts[1];
        if (!xitol.hasOwnProperty(node_name))
            xitol[node_name] = {};
        xitol[node_name][var_type] = {}; //range/significant
        xitol[node_name][var_type].hex_color = parts[2];
        xitol[node_name][var_type].notes = parts[3];
    }
    // console.log(xitol);
    return xitol;
};


parseCSV = function (lines) {
    var palette = ["#FF0000","#00FF00","#0000FF","#FF33CC","#996633","#00CC99","#CC6699","#009933"]
    var table_hist = {"histograms":{}, 'labels':[], 'colors':[]};
    var xitol = {};
    var csv_line = {};
    var rel_freqs_indices = [];
    var line_0 = lines[0];
    var parts_0 = line_0.split("\t");
    for (var idx=0; idx < parts_0.length; idx++)
    {
        if(parts_0[idx] == "By Group Relative Frequency")
            rel_freqs_indices.push(idx)
    }
    var line_l = lines[1];
    var parts_l = line_l.split("\t");
    var color_values = new Set();
    //var line_c = lines[2];
    //var parts_c = line_c.split("\t");
    
    for (var idx = 0; idx < rel_freqs_indices.length; idx++)
    {
        table_hist.labels.push(parts_l[rel_freqs_indices[idx]]);
        table_hist.colors.push(palette[idx%palette.length]);
    }

    for (var line = 3; line <  lines.length; line++)
    {
        var l = lines[line];
        var parts = l.split("\t");
        var node_name = parts[0];
        var is_leaf = parts[1];
        var taxonomy = parts[2];
        var itig = {
            'nats': parseFloat(parts[3]),
            'turnover': parseFloat(parts[4]),
            'pvalue': parseFloat(parts[5]),
            'multtest': parts[6],
            'color': parts[7]
        };
        color_values.add(itig['color']);
        var last_rel_f = rel_freqs_indices[rel_freqs_indices.length-1];
        var itigs = {
            'nats': parseFloat(parts[last_rel_f+1]),
            'turnover': parseFloat(parts[last_rel_f+2]),
            'pvalue': parseFloat(parts[last_rel_f+3]),
            'multtest': parts[last_rel_f+4]
        }
        if (node_name.length>0)
        {
            if (!xitol.hasOwnProperty(node_name))
                xitol[node_name] = {};
            xitol[node_name]["range"] = {"hex_color": itig.color}; //range/clade
            var clade_color = "#000000";
            if(itig.multtest == "True")
                clade_color = "#0000FF";
            xitol[node_name]["clade"] = {"hex_color": clade_color};
            xitol[node_name].taxonomy = taxonomy;
            xitol[node_name].itig = itig;
            xitol[node_name].itigs = itigs;
            xitol[node_name].csv_line = parts.join(";");
            table_hist.histograms[node_name] = [];
            for (var idx = 0; idx < rel_freqs_indices.length; idx++)
                table_hist.histograms[node_name].push(parseFloat(parts[rel_freqs_indices[idx]]));    
        }
    }
    var color_dd = {};
    var colors_weights_dd = {};
    var color_map = Array.from(color_values);
    for (var idx = 0; idx < color_map.length; idx++)
        color_dd[color_map[idx]] = [];
    for (var key in xitol) {
      if (xitol.hasOwnProperty(key)) {
        var cv = xitol[key]["range"]["hex_color"];
        color_dd[cv].push(xitol[key]["itig"]["turnover"]);
      }
    }
    for (var idx = 0; idx < color_map.length; idx++)
    {
        var sum = 0;
        for( var i = 0; i < color_dd[color_map[idx]].length; i++ ){
            sum += color_dd[color_map[idx]][i]; 
        }
        var avg = sum/color_dd[color_map[idx]].length;
        colors_weights_dd[color_map[idx]] = avg;
    }
    var items = Object.keys(colors_weights_dd).map(function(key) {
        return [key, colors_weights_dd[key]];
    });
    items.sort(function(first, second) {
        return first[1] - second[1];
    });
    // array di oggetti da ordinare
    table_hist.color_map = {};
    var temp_lines = [];
    for (var idx = 0; idx < 3; idx++)
        temp_lines.push(lines[idx].split("\t").join(";"));
    table_hist.csv_header = temp_lines.join("\n"); 
    var cnt = 1;
    for (var idx = 0; idx < items.length; idx++)
        if (items[idx][1]!=NaN)
        {
            table_hist.color_map[items[idx][0]] = {};
            table_hist.color_map[items[idx][0]] = {"radius": cnt*3};
            cnt++;
        }
    //table_hist = normalize_hist(table_hist);
    return {'table_hist': table_hist, 'xitol': xitol};
}

function updateLegend(th)
{
    var legend_div = document.getElementById('legend');
    for (var idx = 0; idx < th.colors.length; idx++)
    {
        var color_div = document.createElement('div');
        color_div.style.width = "40px";
        color_div.style.height = "20px";
        color_div.style.float = 'left';
        color_div.style.backgroundColor = th.colors[idx];
        var label_div = document.createElement('div');
        label_div.style.width = "60px";
        label_div.style.float = 'left';
        label_div.innerHTML = th.labels[idx];
        var brel = document.createElement('br');
        brel.style.clear = 'left';
        legend_div.appendChild(color_div);
        legend_div.appendChild(label_div);
        legend_div.appendChild(brel);
    }
}

function parseCSVLines(lines)
{
    var dd = parseCSV(lines);
    opts.table_hist = dd.table_hist;
    opts.xitol = dd.xitol;
    document.getElementById("csv-area").value = dd.table_hist.csv_header + "\n"; 
    //updateTree(tree,opts);
    updateLegend(opts.table_hist);
}


window.onload = function () {
    var fileInputCSV = document.getElementById('fileInputCSV');
    var fileInputTree = document.getElementById('fileInputTree');
    fileInputCSV.addEventListener('change', function(e) {
        var num_files = fileInputCSV.files.length;
        if ( (num_files>0) && (fileInputCSV.files[0].name.indexOf('.mibybranch')>-1) )
        {
            var file = fileInputCSV.files[0];
            var reader = new FileReader();
            reader.onload = function (e)
            {
                console.log(reader.result);
                var lines = reader.result.split('\n');
                //debugger;
                parseCSVLines(lines);
                fileInputTree.disabled = false;
            }
            reader.readAsText(file);
        }
    });
    fileInputTree.addEventListener('change', function(e) {
        var num_files = fileInputTree.files.length;
        if ( (num_files>0) && (fileInputTree.files[0].name.indexOf('.TreeLabeled')>-1) )
        {
            var file = fileInputTree.files[0];
            var reader = new FileReader();
            reader.onload = function (e)
            {
                var lines = reader.result.split('\n');
                opts.tree.original = lines[0];
                opts.tree.data = lines[0];
                var start = new Date().getTime();
                updateTree(tree,opts);
                /*var mynode = d3.select(".root");
                mynode.append("circle")
                        .attr("x",1000)
                        .attr("y",1000)
                        .attr("r",800+30)
                        .attr("stroke", "black")
                        .attr("stroke-width",0.5)
                        .attr('fill', 'none');
                */
                var end = new Date().getTime();
                var time = end - start;
                console.log('Execution time load new tree: ' + time);
                //debugger;
                computeInnerLeavesHistograms(tree.data());
            }
            reader.readAsText(file);
        }
    });
};

/*
if (file.name.indexOf('.TreeLabeled')==-1)
                    {
                        if (lines[0].split(',')[0]=="LABELS")
                        {
                            //histxitol
                            opts.table_hist = parseTableHist(lines);
                            updateTree(tree,opts);
                        }
                        else
                        {
                            opts.xitol = parseXItol(lines);
                            updateTree(tree, opts);
                        }
                    }
                    else
                    {
                        // tree in newick format
                        opts.tree.original = lines[0];
                        opts.tree.data = lines[0];
                        updateTree(tree,opts); 
                        debugger;
                        computeInnerLeavesHistograms(tree.data());   
                    }    
*/
// update the elements
function updateCircleRadius(nRadius) {
    // adjust the text on the range slider
    d3.select("#nRadius-value").text(nRadius);
    d3.select("#nRadius").property("value", nRadius);
    // update the circle radius
    d3.selectAll("circle") 
        .attr("r", nRadius);
}

function initHoverTooltip() {
    var div = d3.select("body").append("div") 
        .attr("id","tooltip")  
        .attr("class", "tooltip")               
        .style("opacity", 0);    
    var myElements = document.querySelectorAll(".tooltip");
    for (var i = 0; i < myElements.length; i++) {
        var st = myElements[i].style;
        st.position = "absolute";
        st.textAlign = "center";
        st.width = "100px";
        st.height = "28px";
        st.padding ="2px";
        st.font = "12px sans-serif"; 
        st.background = "lightsteelblue";
        st.border = "0px";
        st.borderRadius = "8px";
        st.pointerEvents = "none";
    }
}

function testRootChain() {
    var dataset = [
    {
        "name": "root",
        "id":1
    },
    {
        "name": "child1",
        "id":2
    },
    {
        "name": "child1.1",
        "id":3
    }];
    
}

function testContextMenu()
{
    var menu = [
            {
                title: 'Item #1',
                action: function(elm, d, i) {
                    console.log('Item #1 clicked!');
                    console.log('The data for this circle is: ' + d);
                }
            },
            {
                title: 'Item #2',
                action: function(elm, d, i) {
                    console.log('You have clicked the second item!');
                    console.log('The data for this circle is: ' + d);
                }
            }
        ]
        var data = [1, 2, 3];
        var g = d3.select('#context-menu').append('svg')
            .attr('width', 200)
            .attr('height', 400)
            .append('g');
        g.selectAll('circles')
            .data(data)
            .enter()
            .append('circle')
            .attr('r', 30)
            .attr('fill', 'steelblue')
            .attr('cx', function(d) {
                return 100;
            })
            .attr('cy', function(d) {
                return d * 100;
            })
            .on('contextmenu', d3.contextMenu(menu, {
                onOpen: function() {
                    console.log('opened!');
                },
                onClose: function() {
                    console.log('closed!');
                }
            })); // attach menu to element
}