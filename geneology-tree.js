
function Geneology(pop) {
    var _size;
    function size(size) {
	if (!arguments.length) return _size;
	_size = size;
	return this;
    };
    function nodes() {
	// flatten the geneology into a single array, storing the generation
	// with the individual
	var flat = [], tmp;
	for (var gen = 1; gen < pop.gens(); gen++) {
	    for (var i = 0; i < pop.individuals[gen].length; i++) {
		ind = pop.individuals[gen][i];
		var obj = {gen: gen,
			   loci: ind.loci,
			   id: ind.id,
			   mid: ind.mid,
			   pid: ind.pid,
			   children: ind.children};
		flat.push(obj);
	    }
	}
	return flat;
    };
    function edges() {
	var edges = [];
	for (var gen = 2; gen < pop.gens(); gen++) {
	    // skip first generation (no parents)
	    pop.individuals[gen].forEach(function(ind) {
		ind.children.forEach(function(kid) {
		    edges.push({gen:gen, id:ind.id, child: kid});
		});
	    });
	}
	return edges;
    };
    return {pop: pop,
	    size: size,
	    edges: edges,
	    nodes: nodes};
	    
}

// example sim
pop = Population(100, constantLinkage(100, 0.01));

dem = Demography().popSizeChangeEvent(15, 10).popSizeChangeEvent(20, 10);
//dem = Demography().popSizeChangeEvent(7, 10);

// side effects: sim (may change) TODO
DiploidWrightFisher(pop, dem); // runs pop through demography

var g = Geneology(pop);

// viz
var margin = {top: 40, right: 120, bottom: 20, left: 120},
    width = 2000 - margin.right - margin.left,
    height = 1900 - margin.top - margin.bottom;

var svg = d3.select("body").append("svg")
	.attr("width", width + margin.right + margin.left)
	.attr("height", height + margin.top + margin.bottom);

// build the tree
var nodes = svg.selectAll("circle");

// TODO write node to position function

function indXPos(id) {
    return id * 30 + 10; 
};

function indYPos(gen) {
    return gen * 50 + 10; 
};

var node_data = g.nodes();
// node_data = [
//     {id: 1, mid: 1, pid: 4, children: [4]}
// ];

function alleleFiller(svg, col1, col2) {
    // append allele color scheme to svg, return func to make colors
    var defs = svg.append("defs");
    var grad_01 = defs.append("linearGradient").attr("id", "allele_01");
    grad_01.append("stop").attr("offset", "0%").style("stop-color", col1);
    grad_01.append("stop").attr("offset", "50%").style("stop-color", col1);
    grad_01.append("stop").attr("offset", "50%").style("stop-color", col2);
    grad_01.append("stop").attr("offset", "100%").style("stop-color", col2);
    var grad_00 = defs.append("linearGradient").attr("id", "allele_00");
    grad_00.append("stop").style("stop-color", col1);
    var grad_11 = defs.append("linearGradient").attr("id", "allele_11");
    grad_11.append("stop").style("stop-color", col2);
    var fill = ["allele_00", "allele_01", "allele_11"];
    fill = fill.map(function(x) { return "fill:url(#" + x + ")"; });
    function fill_func(mallele, pallele) {
	return fill[mallele + pallele];
    };
    return fill_func;
}

var alleleFill = alleleFiller(svg, "#d7191c", "#2c7bb6");

nodes.data(node_data).enter().append("circle")
    .attr("cx", function(x) { return indXPos(x.id); })
    .attr("cy", function(x) { return indYPos(x.gen); }) 
    .attr("r", "5")
    //.attr("style", "fill:steelblue;");
    .attr("style", function(x) {
	return alleleFill(x.loci[0][0], x.loci[1][0]) + ";";
    });

var edge_data = g.edges();
// edge_data = [];

nodes.data(edge_data).enter().append("line")
    .attr({
	x1: function(n) { return indXPos(n.id); },
	x2: function(n) { return indXPos(n.child); },
	y1: function(n) { return indYPos(n.gen-1); },
	y2: function(n) { return indYPos(n.gen); }
    })
    .style("fill", "none")
    .style("stroke", "gray")
    .style("stroke-opacity", 0.6);
	





