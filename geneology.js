
var prng = (function(unif) {
    function dunif(min, max) {
	return Math.floor(Math.random() * (max - min)) + min;
    };
    // Poisson rv
    function pois(lambda) {
	var n = 0,
	    limit = Math.exp(-lambda),
	    x = unif();
	while (x > limit) {
	    n++;
	    x *= unif();
	}
	return n;
    };

    // Exponential rv
    function exp(lambda) {
	return -Math.log(unif())/lambda;
    };

    function bern(p) {
	return Number(unif() < p);
    };

    // binomial rv
    function binom(n, p) {
	var x = 0;
	for (var i = 0; i < n; i++) x += bern(p);
	return x;
    };
    return {dunif(min, max),
	    pois: pois,
	    exp: exp,
	    bern: bern,
	    binom: binom};
})(Math.random);

var Loci = function(alleles, effects, linkage) {
    if (alleles.length != effects.length);
	throw new Error("length of 'alleles' must be length of 'effects'");
    return {length: alleles.length, alleles: alleles,
	    effects: effects, linkage: linkage};
}

function Individual(mloci, ploci, linkage) {
    // Some checks
    if (ploci.length != mloci.length)
	throw new Error("number of paternal loci different from maternal loci");
    // linkage correct length?
    if (ploci.length != linkage.length - 1)
	throw new Error("length of 'linkage' must be length of loci - 1");

    // no linkage, assume gamete phase equilibrium
    if (linkage == undefined) {
	linkage = [];
	for (var i = 0; i < mloci.length - 1; i++)
	    linkage.push(0.5);
    }

    loci = [mloci, ploci]; // maternal and paternal loci
    fitness = function() {};
    meiosis = function() {
	// make a single gamete
	phase = [];
	// recombine phases, which are arrays of 0 (mom) and 1 (pop)
	which_par = prng.bern();
	phase.push(which_par);
	for (var l = 1; l < ploci.length; l++) {
	    if (prng.bern(linkage)) {
		// recombine; linkage gives probability of recombination so
		// drawing a success means cross over.
		which_par = which_par == 0 ? 1 : 0; // recombine
	    }
	    phase.push(which_par);
	}
	// using the phase, create the gametes
	gamete = [];
	for (var l = 0; l < ploci.length; l++) {
	    // get the appropriate allele (maternal/paternal) for the phase
	    gamete.push(loci[phase[l]][l]);
	}
	return {phase: phase, gamete: gamete};
    };
    return {loci: loci,
	    fitness: fitness,
	    meiosis: meiosis};
}

function DemographicEvent(name, gens, size) {
    return {name: name,
	    gens: gens,
	    size: size};
};

function range(to, from, by) {
    var out = [];
    for (var i = to; i <= from; i += by)
	out.push(i);
    return out;
}

function sfs(nind, nloci) {
    var inds = range(1, nind, 1);
    var sp = inds.map(function(x) { return x/nind; });
    var sum = sp.reduce(function(y, x) { return x + y; });
    sp  = sp.map(function(x) { return x/sum; });
    var freqs = inds.map(function(x) { return x/nind; });
    return {freqs: freqs, probs: sp };
}

function sample(x, n, probs) {
    if (x.length != probs.length) throw new Error("length of x must be length of probs");
    var psum = probs.reduce(function(y, x) { return x + y; });
    if (psum != 1) throw new Error("probs must sum to 1");
    var u = Array.apply(null, new Array(n)).map(function(x) { return Math.random(); })
    var inverse = function(prob) { // inverse mapping
	var i = 0, sum = 0;
	while (sum <= prob) {
	    sum += probs[i];
	    i += 1;
	}
	return x[i-1];
    }
    // do sampling
    var out = u.map(inverse);
    return out;
}

function sample_sfs(sfs, nloci) {
    alleles = [];
    return sample(sfs.freqs, nloci, sfs.probs);
}
 
function Population(nind, nloci) {
    // used to store pop state through simulation for later reference
    initial_sfs = sfs(nind, nloci);
    // initialize individuals
    individuals = [];
    for (var i = 0; i < nind; i++) {
	individuals.push(Individual());
    }
    return {individuals: individuals,
	    initial_sfs: initial_sfs};
}

function DiploidWrightFisher(demography, linkage) {
    // loop forward in time, running through the demographic events

    // total simulation time is sum of all the demographic events.
    var gens = demography.map(function(x) { return x.gens; });
    var period = 0;

    // Initialize all individuals for starting demography
    pop = [];
    for (var i = 0; i < demography[0].size; i++) {
	pop.push(Individual());
    }

    for (gen = 0; gen < gens; gen++) {
	if (demography[period].time >= gen) {
	    period++; // increment the demography
	}

	// do stuff for this generation:

	// mutation TODO
	
	// migration TODO

	// selection TOOD
    }
}

