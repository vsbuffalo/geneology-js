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

