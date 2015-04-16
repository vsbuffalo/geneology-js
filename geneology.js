EPS = 0.000001
// simple sampling distributions

var Dist = (function(unif) {
    // discrete uniform
    function dunif(min, max) {
	return Math.floor(unif() * (max - min)) + min;
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
    return {dunif: dunif,
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

function isValidLinkage(loci, linkage) {
    // linkage correct length?
    if (linkage == undefined || loci == undefined ||
	loci.length - 1 != linkage.length)
	throw new Error("length of 'linkage' must be length of loci - 1 (loci: "
			+ loci.length + "; linkage: " + linkage.length + ")");
}

function Individual(mloci, ploci, id, mid, pid) {
    // Some checks
    if (ploci.length != mloci.length)
	throw new Error("number of paternal loci different from maternal loci");
    loci = [mloci, ploci]; // maternal and paternal loci
    fitness = function() {};
    meiosis = function(linkage) {
	isValidLinkage(ploci, linkage);
	// make a single gamete
	phase = [];
	// recombine phases, which are arrays of 0 (mom) and 1 (pop)
	which_par = Dist.bern();
	phase.push(which_par);
	for (var l = 0; l < ploci.length; l++) {
	    if (Dist.bern(linkage[l])) {
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
	    meiosis: meiosis,
	    id: id,
	    mid: mid,
	    pid: pid};
}

function DemographicEvent(name, gens, size) {
    return {name: name,
	    gens: gens,
	    size: size};
};

function range(to, from, by) {
    var out = [];
    for (var i = to; i <= from; i += by) out.push(i);
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
    if (Math.abs(psum - 1) > EPS) throw new Error("probs must sum to 1 +/- " + EPS + " (sum " + psum + " )");
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

function constantLinkage(nloci, r) {
    // no linkage, assume gamete phase equilibrium
    return Array.apply(null, new Array(nloci-1))
	.map(function(x) { return r; })
}

function sample_sfs(sfs, nloci) {
    alleles = [];
    return sample(sfs.freqs, nloci, sfs.probs).
	map(function(freq) {return Number(Dist.bern(freq));});
}
 
function Population(nloci, linkage) {
    // used to store pop state through simulation for later reference
    // Should contain subpopulations/demes
    var individuals = [];
    var demes = {}; // TODO
    var history = [];
    var initial_sfs;
    function init(nind) {
	if (nind == undefined || nind < 2)
	    throw new Error("nind must be > 2");
	initial_sfs = sfs(nind, nloci); // TODO resolution? multiply nind by 1000?
	// initialize individuals
	var first_gen = [], mom, pop;
	for (var i = 0; i < nind; i++) {
	    mom = sample_sfs(initial_sfs, nloci);
	    pop = sample_sfs(initial_sfs, nloci);
	    first_gen.push(Individual(mom, pop, i, null, null));
	}
	individuals.push(first_gen);
	return this;
    }
    function last_gen() {
	// return the last generation
	return individuals[individuals.length-1];
    }
    function random_individual() {
	// sample a random individual from the last generation
	size = individuals[individuals.length-1].length;
	return individuals[individuals.length-1][Dist.dunif(0, size-1)];
    }
    function mate(popsize) {
	// make the next generation
	var new_gen = [];
	for (var i = 0; i < popsize; i++) {
	    var mom = random_individual();
	    var pop = random_individual();
	    var kid = Individual(mom.meiosis(linkage).gamete,
				 pop.meiosis(linkage).gamete,
				 i, mom.id, pop.id);
	    new_gen.push(kid);
	}
	if (new_gen.length > 0)
	    individuals.push(new_gen);
	return this;
    }
 
    function gens() { return individuals.length; };

    return {init: init,
	    random_individual: random_individual,
	    mate: mate,
	    gens: gens,
	    last_gen: last_gen,
	    individuals: individuals,
	    initial_sfs: initial_sfs};
}

function Demography() {
    var events = [];
    function popSizeChangeEvent(size, gens, name) {
	events.push({size: size, gens: gens, name: name});
	return this;
    }
    return {events: events, popSizeChangeEvent: popSizeChangeEvent};
}

function DiploidWrightFisher(pop, demography) {
    // Push all populations through demography, each pop is a simulation.
    dem_events = demography.events;

    // total simulation time is sum of all the demographic events.
    var gens = dem_events.map(function(x) { return x.gens; });
    var period = 0;

    // get first population size initialize population
    var init_nind = dem_events[0].size;
    pop.init(init_nind);

    // Initialize all individuals for starting demography
    for (gen = 0; gen < gens; gen++) {
	if (dem_events[period].time >= gen) {
	    period++; // increment the demography
	}
	// get population size for tihs generation
	popsize = dem_events[period].size;
	pop.mate(popsize);

	// do other stuff for this generation:

	// mutation TODO
	
	// migration TODO

	// selection TOOD
    }
}



pop = Population(100, constantLinkage(100, 0.01))

dem = Demography().popSizeChangeEvent(15, 10)

wf = DiploidWrightFisher(pop, dem)
