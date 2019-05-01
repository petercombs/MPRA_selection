var nodeStyler = function (element, data) {
	if(data.hasOwnProperty('internal_label')) {
		// If the node has an internal label (see below), we display it
		// next to the node by creating a new 'text' element.

		// Make sure we don't already have an internal label node on this SVG node!
		var label = element.selectAll(".internal_label");
		if(label.empty()) {
			var text_label = element.append("text");

			// TODO: Once we're happy with how these elements look,
			// we should move all this complexity into CSS classes.
			text_label.classed("internal_label", true)
				.text(data.internal_label)
				.attr("dx", ".4em")
				.attr("dy", ".3em")
				.style("font-style", "italic")
				.attr("text-anchor", "start")
				.attr("alignment-baseline", "middle");

			// If the internal label has the same label as the currently
			// selected phyloreference, make it bolder and turn it blue.
			/*
			if(
				selected_phyloref !== undefined &&
				selected_phyloref.hasOwnProperty('label') &&
				selected_phyloref.label == data.internal_label
			) {
				text_label.style('fill', 'blue')
					.style('font-weight', 'bolder');
			}
			*/
		}
	}
}


//example_tree = " (((((Papio_anubis:0.005152,Macaca_fascicularis:0.011288)N5:0.063148,(((Gorilla_gorilla:0.005757,(Homo_sapiens:0.002882,(Pan_paniscus:0.003623,Pan_troglodytes:0.000000)N10:0.000000)N9:0.007273)N8:0.000000,Pongo_abelii:0.017407)N7:0.000000,Nomascus_leucogenys:0.014595)N6:0.010172)N4:0.000489,(Callithrix_jacchus:0.008288,Saimiri_boliviensis:0.025922)N11:0.029016)N3:0.099973,Otolemur_garnettii:0.220235)N2:0.052110,Lipotes_vexillifer:0.146237)N1;";
example_tree = "";
d3.text('ECR11-smith.tree', function(d){
	tree(example_tree=d.replace(/:[0-9.]*/gi, ':1'))
		.options(
			{
				"reroot":false,
				//"draw-size-bubbles": true,
				"node-circle-size": (d)=>(2),
			})
		.spacing_x(30)
		.spacing_y(30)
		.layout();
	;
	_.each(tree.get_nodes(), function(node) {
		if(node.children && node.name.startsWith("N")) {
			node.internal_label = node.name;
		}
	});

	tree
		.placenodes()
		.update()
	;
	tree.selection_callback(update_selection);
});

example_tree = example_tree.replace(/:[0-9.]*/gi, ':1');
var tree = d3.layout.phylotree()
	.svg(d3.select("#tree_display"))
	.style_nodes(nodeStyler)
;

//tree(example_tree)
//tree


all_seqs = {};
current_seqs = {};

sel = [];
dataset = [];

var colors = {
"A" : "red",
	"G" : "yellow",
	"C" : "blue",
	"T" : "cadetblue",
	"N" : "gray",
	"-" : "lightgray"
};
gridSize = 15;
pad =  10*gridSize;

alignment = d3.select("#alignment_display");

sampleLabels = [];
points = [];

d3.json("ECR11-smith.json", function(error, data){
	all_seqs = data;
	current_seqs = all_seqs;
	points = alignment.selectAll(".point");
	points.data([]).exit().remove();
	sampleLabels = alignment.selectAll(".sampleLabel");
	sampleLabels.data([]).exit().remove();
	sampleLabels.data([]).enter().append("text");
	update_selection([])


}
);

var tooltip = d3.select("body").append("div")
		.attr("class", "tooltip")
		.style("opacity", 0)
	;

function update_selection(seqnames){
	seqnames = seqnames.map(d=>d.name);
	seqnames2 = []
	console.log(seqnames);
	seqnames.forEach(function(name){
		console.log(name);
		console.log(tree.get_node_by_name);
		parname = tree.get_node_by_name(name).parent.name;
		if (parname in seqnames){
			seqnames2.push(name);
		}else{
			seqnames2.push(parname);
			seqnames2.push(name);
		}
	});
	seqnames = seqnames2;

	if (seqnames.length == 0)
	{
		current_seqs = all_seqs;
	}else{
		seqnames.push("original")
		current_seqs = _.pick(all_seqs, seqnames);
	}
	sel = current_seqs;
	[dataset, variable_pos] = collapse_sequences(current_seqs);
	breaks = []
	for (i = 0; i < variable_pos.length-1; i++){
		if (variable_pos[i+1] - variable_pos[i] > 1)
			breaks.push(i)
	}

	//alignment.selectAll('.breakline').data([]).exit().remove()
	breaklines = alignment.selectAll(".breakline").data(breaks, d=>d);
	breaklines.enter().append('line')
		.attr('class', 'breakline')
		.attr('x1', d=> pad + (d+1) * gridSize)
		.attr('x2', d=>pad + (d+1) * gridSize)
		.attr('y1', 0)
		.attr('y2', gridSize*(_.keys(current_seqs).length+1))
		.style('stroke', 'black')
		.style('stroke-width', 2);

	breaklines.transition()
		.duration(200)
		.attr('x1', d=> pad + (d+1) * gridSize)
		.attr('x2', d=>pad + (d+1) * gridSize)
		.attr('y1', 0)
		.attr('y2', gridSize*(_.keys(current_seqs).length+1))
	;

	breaklines.exit().remove();


	//points.data([]).exit().remove();
	points = alignment.selectAll(".point").data(dataset, d=>d.seqname + d.pos);

	points.transition()
			.duration(500)
			.attr("x", function(d){ return pad + d.row * gridSize+1})
			.attr("y", function(d){ return d.col * gridSize+1})
			.style("fill", function(d){ return colors[d.base]})
	;

	points
		.enter()
		.append("rect")
			.attr("class", "rect bordered point")
			.attr("x", function(d){ return pad + d.row * gridSize+1})
			.attr("y", function(d){ return d.col * gridSize+1})
			.attr("width", gridSize-2)
			.attr("height", gridSize-2)
			.style("fill", function(d){ return colors[d.base]})
			.on("mouseover", function(d){
				tooltip.transition()
					.duration(200)
					.style("opacity", .9);
				tooltip.html(  "col - " +
					(d.pos + 1) + "<br/>" +  "base - " + d.base + '<br /> seq - ' + d.seqname + "<br /> seqpos - " + d.seqpos)
					.style("left", (d3.event.pageX + 5) + "px")
					.style("top", (d3.event.pageY - 10) + "px");
			})
			.on("mouseout", function(d){
				tooltip.transition()
					.duration(500)
					.style("opacity", 0)
			});


	points.exit().remove();

	//sampleLabels.data([]).exit().remove();

	sampleLabels = alignment.selectAll(".sampleLabel");
	sampleLabels = sampleLabels.data(_.keys(current_seqs), d=>d);
	sampleEnter = sampleLabels
		.enter().append("text")
		.attr('class', 'sampleLabel')
		.attr("x",pad)
		.attr("y", function(d,i){ return i * gridSize;})
		.style("text-anchor", "end")
		.attr("transform", "translate(-6," + gridSize + ")" )
		.text(d=>d)
	;
	sampleLabels
		.transition()
		.duration(200)
		.attr("y", function(d,i){ return i * gridSize;})
		.text(function(d,i){ return (d); })
	;
	sampleLabels.exit().remove()

}

function collapse_sequences(seq_obj){
	seqnames = _.keys(seq_obj);
	sel = seqnames;
	seqlen = seq_obj[seqnames[0]].length;
	variable_pos = [];
	for (pos = 0; pos < seqlen; pos++){
		bases = new Set();
		seqnames.forEach(function(d, i){
			if (d != 'original')
				bases.add(seq_obj[d].charAt(pos));
		});

		bases.delete('N');
		// Ignore masked bases
		if (bases.size > 1){
			variable_pos.push(pos);
		}
	}

	var dataset = [];
	for (var i=0; i<seqnames.length; i++){
		var key = seqnames[i];
		var x = seq_obj[key];
		for (var j=0; j<variable_pos.length; j++){
			pos = variable_pos[j];
			orig_pos = x.substring(0, pos).replace(/-/g, '').length
			var entry = {"base" : x[pos], "row" : j, "col" : i, "seqname": key, "pos": pos, "seqpos": orig_pos+1};
			dataset.push(entry);
		};
	};
	return [dataset, variable_pos];

}



