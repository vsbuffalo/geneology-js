EPS = 0.000001;
// simple sampling distributions

// TODO: Individual shouldn't have loci/phase unpacked in constructor

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

function Individual(mloci, ploci, id, mid, pid, mphase, pphase) {
    // Some checks
    if (ploci.length != mloci.length)
	throw new Error("number of paternal loci different from maternal loci");
    var loci = [mloci, ploci]; // maternal and paternal loci
    var phases = [mphase, pphase];
    var children = [], is_mom = [];
    function fitness() {};
    function meiosis(linkage) {
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
    function addChild(kid, is_mom) {
	children.push(kid.id);
	this.is_mom.push(Number(is_mom));
	return this;
    };
    return {loci: loci,
	    fitness: fitness,
	    meiosis: meiosis,
	    id: id,
	    mid: mid,
	    pid: pid,
	    phases: phases,
	    children: children,
	    is_mom: is_mom,
	    addChild: addChild};
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
    var u = Array.apply(null, new Array(n)).map(function(x) { return Math.random(); });
    var inverse = function(prob) { // inverse mapping
	var i = 0, sum = 0;
	while (sum <= prob) {
	    sum += probs[i];
	    i += 1;
	}
	return x[i-1];
    };
    // do sampling
    var out = u.map(inverse);
    return out;
}

function constantLinkage(nloci, r) {
    // no linkage, assume gamete phase equilibrium
    return Array.apply(null, new Array(nloci-1))
	.map(function(x) { return r; });
}

function sample_sfs(sfs, nloci) {
    alleles = [];
    return sample(sfs.freqs, nloci, sfs.probs).
	map(function(freq) {return Number(Dist.bern(freq));});
}
 
function Population(nloci, linkage, mutate) {
    // used to store pop state through simulation for later reference
    // Should contain subpopulations/demes
    var individuals = [];
    var demes = {}; // TODO
    var history = [];
    var initial_sfs;
    var mutate = mutate;
    function init(nind) {
	if (nind == undefined || nind < 2)
	    throw new Error("nind must be > 2");
	initial_sfs = sfs(nind, nloci); // TODO resolution? multiply nind by 1000?
	// initialize individuals
	// first invididual 'eve' is just place holder, the parent of all first gens
	// makes drawing trees easier (since there's a root).
	var first_gen = [], mom, dad, kid, mom_phase, dad_phase;
	for (var i = 0; i < nind; i++) {
	    mom = sample_sfs(initial_sfs, nloci);
	    dad = sample_sfs(initial_sfs, nloci);
	    mom_phase = Array.apply(null, new Array(nloci)).map(function(){return 1;});
	    dad_phase = Array.apply(null, new Array(nloci)).map(function(){return 0;});
	    kid = Individual(mom, dad, i, 0, 0, mom_phase, dad_phase);
	    first_gen.push(kid);
	}
	individuals.push(first_gen);
	return this;
    }
    function last_gen() {
	// return the last generation
	return individuals[individuals.length-1];
    }
    function random_individual(except) {
	// sample a random individual from the last generation
	// ignore is an individual to ignore (e.g. the current individual)
	var max = individuals[individuals.length-1].length - 1;
	var sample_space = individuals[individuals.length-1];
	if (typeof except != undefined) {
	    sample_space = sample_space.filter(function(x) {
		return x.id != except;
	    });
	    max--;
	}
	return sample_space[Dist.dunif(0, max)];
    }
    function mate(popsize, selfing) {
	// make the next generation
	var new_gen = [];
	for (var i = 0; i < popsize; i++) {
	    var mom = random_individual();
	    var except = typeof selfing == undefined || !selfing ? mom.id : undefined;
	    var dad = random_individual(mom.id);
	    var mom_gamete = mom.meiosis(linkage),
		dad_gamete = dad.meiosis(linkage);
	    mom_gamete.gamete = mutate(mom_gamete.gamete);
	    dad_gamete.gamete = mutate(dad_gamete.gamete);
	    var kid = Individual(mom_gamete.gamete,
				 dad_gamete.gamete,
				 i, mom.id, dad.id,
				 mom_gamete.phase,
				 dad_gamete.phase);
	    mom.addChild(kid, true);
	    dad.addChild(kid, false);
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
	    mutate: mutate,
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

function bernMutater(prob) {
    function mutater(gamete) {
	return gamete.map(function(a) {
	    var mutate = Dist.bern(prob);
	    if (!mutate) return a;
	    return a == 0 ? 1 : 0;
	});
    };
    return mutater;
}

function DiploidWrightFisher(pop, demography) {
    // Push all populations through demography, each pop is a simulation.
    dem_events = demography.events;

    // total simulation time is sum of all the demographic events.
    var gens = d3.sum(dem_events.map(function(x) { return x.gens; }));
    var max_size = Math.max.apply(null, dem_events.map(function(x) {return x.size;}));
    var period = 0;

    // get first population size initialize population
    var init_nind = dem_events[0].size;
    pop.init(init_nind);

    // Initialize all individuals for starting demography
    var event_gen = 0;
    for (gen = 0; gen < gens; gen++) {
	if (dem_events[period].gens <= event_gen) {
	    period++; // increment the demography
	    event_gen = 0;
	}
	// get population size for tihs generation
	popsize = dem_events[period].size;
	pop.mate(popsize);

	// do other stuff for this generation:
	
	// migration TODO

	// selection TOOD
	event_gen += 1; // which generation we are in demographic event
    }
    return {pop: pop, gens: gens, max_size: max_size};
}

