// settings
gen_x_buffer = 40; // multiplier (determines node horizontal spacing)
gen_y_buffer = 40; // multiplier (determines node vertical spacing)
gen_x_margin = 10; // margin (horizontal margin between nodes)
gen_y_margin = 10; // margin (vertical margin between nodes)
allele_x_buff = 11;
allele_colors = ["#d8b365", "#5ab4ac"];

function Geneology(pop) {
    var _size;
    function size(size) {
	if (!arguments.length) return _size;
	_size = size;
	return this;
    };
    function nodes() {
	// nodes based on individuals
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
    function backtrace(ind) {

    };
    function forwardtrace(ind) {

    };
    function edges() {
	var edges = [];
	for (var gen = 1; gen < pop.gens(); gen++) {
	    // skip first generation (no parents)
	    pop.individuals[gen].forEach(function(ind) {
		ind.children.forEach(function(kid, kid_i) {
		    // TODO: beware fencepost error here?
		    var is_mom = ind.is_mom[kid_i];
		    var kid_phase = pop.individuals[gen][ind.id].phases[is_mom];
		    edges.push({gen:gen+1,
				id:ind.id,
				child: kid,
				is_mom: is_mom,
				child_phase: kid_phase});
		});
	    });
	}
	return edges;
    };
    return {pop: pop,
	    size: size,
	    edges: edges,
	    nodes: nodes,
	    forwardtrace: forwardtrace,
	    backtrace: backtrace};
	    
}

function svgSimSetup(max_size, gens) {
    // take a wright fisher sim and produce width/height
    return {
	width: indXPos(max_size) + gen_x_margin*2,
	height: indYPos(gens) + gen_y_margin*2
    };
}

// example sim
var mutater = bernMutater(0.01);
var pop = Population(100, constantLinkage(100, 0.01), mutater);

var dem = Demography().popSizeChangeEvent(20, 10)
	.popSizeChangeEvent(4, 10)
	.popSizeChangeEvent(10, 16);
//dem = Demography().popSizeChangeEvent(7, 10);

// side effects: sim (may change) TODO
var wf = DiploidWrightFisher(pop, dem); // runs pop through demography

var g = Geneology(pop);

// viz
var setup = svgSimSetup(wf.max_size, wf.gens);
var margin = {top: 40, right: 120, bottom: 20, left: 120},
    width = setup.width - margin.right - margin.left,
    height = setup.height - margin.top - margin.bottom;

var svg = d3.select("body").append("svg")
	.attr("width", width + margin.right + margin.left)
	.attr("height", height + margin.top + margin.bottom);

// build the tree
var nodes = svg.selectAll("circle");

function indXPos(id) {
    return id * gen_x_buffer + gen_x_margin; 
};

function indYPos(gen) {
    return gen * gen_y_buffer + gen_y_margin; 
};

function alleleXPos(id, haplo) {
    // haplo is an offset, [is_mom, kidphase]
    var offset = haplo == undefined ? 0 : haplo;
    return id * gen_x_buffer + haplo*allele_x_buff + gen_x_margin; 
};

var node_data = g.nodes();
// node_data = [
//     {id: 1, mid: 1, pid: 4, children: [4]}
// ];

function genotypeFiller(svg, col1, col2) {
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

function alleleFiller(svg, col1, col2) {
    var defs = svg.append("defs");
    var grad_0 = defs.append("linearGradient").attr("id", "allele_0");
    grad_0.append("stop").style("stop-color", col1);
    var grad_1 = defs.append("linearGradient").attr("id", "allele_1");
    grad_1.append("stop").style("stop-color", col2);
    var fill = ["allele_0", "allele_1"];
    fill = fill.map(function(x) { return "fill:url(#" + x + ")"; });
    function fill_func(allele) {
	return fill[allele];
    };
    return fill_func;
}

var genotypeFill = genotypeFiller(svg, allele_colors[0], allele_colors[1]);
var alleleFill = alleleFiller(svg, allele_colors[0], allele_colors[1]);

var edge_data = g.edges();
// edge_data = [];

// === Individual nodes
// nodes.data(edge_data).enter().append("line")
//     .attr({
// 	x1: function(n) { return indXPos(n.id); },
// 	x2: function(n) { return indXPos(n.child); },
// 	y1: function(n) { return indYPos(n.gen-1); },
// 	y2: function(n) { return indYPos(n.gen); }
//     })
//     .style("fill", "none")
//     .style("stroke", "gray")
//     .style("stroke-opacity", 0.4);
	
// nodes.data(node_data).enter().append("circle")
//     .attr("id", function(x) { return x.gen + "-" + x.id; })
//     .attr("cx", function(x) { return indXPos(x.id); })
//     .attr("cy", function(x) { return indYPos(x.gen); }) 
//     .attr("r", "5")
//     //.attr("style", "fill:steelblue;");
//     .attr("style", function(x) {
// 	return alleleFill(x.loci[0][0], x.loci[1][0]) + ";" +
// 	    // TODO broken:
// 	    "fill-opacity:"+[0.3, 1][Number(x.children.length > 0)] + ";";
//     });

// === Allele nodes

nodes.data(node_data).enter().append("g")
    .selectAll(".allele")
    .data(function(ind) {
	var locus = 0;
	var obj = [{id: ind.id, gen: ind.gen, is_mom: 1, allele: ind.loci[0][locus]},
 		   {id: ind.id, gen: ind.gen, is_mom: 0, allele: ind.loci[1][locus]}];
	//debugger;
	return obj;
    })
    .enter().append("circle")
    .attr("cx", function(x) { return alleleXPos(x.id, x.is_mom); })
    .attr("cy", function(x) { return indYPos(x.gen); }) 
    .attr("r", "4")
    .attr("style", function(x) {
	return alleleFill(x.allele) + ";";
    });

svg.selectAll(".link")
    .data(edge_data)
    .enter().append("line")
    .attr({
	x1: function(n) {
	    var locus = 0;
	    // TODO: need to invert phase here -- why?
	    return alleleXPos(n.id, Number(!n.child_phase[locus]));
	},
	x2: function(n) {
	    return alleleXPos(n.child, n.is_mom);
	},
	y1: function(n) { return indYPos(n.gen-1); },
	y2: function(n) { return indYPos(n.gen); }
    })
    .style("fill", "none")
    .style("stroke", "gray")
    .style("stroke-opacity", 0.4);
	
