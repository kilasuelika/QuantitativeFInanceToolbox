*Calculate emperical cdf;
proc sort data=grad.logcurr5eth;
   by btc;
run;
data grad.logcurr5eth;
set grad.logcurr5eth nobs=totalobs;
   ecdf_btc = _n_ / totalobs;
run;
proc sort data=grad.logcurr5eth;
   by eth;
run;
data grad.logcurr5eth;
set grad.logcurr5eth nobs=totalobs;
   ecdf_eth = _n_ / totalobs;
run;
proc sort data=grad.logcurr5eth;
	by date;
run;
*Remove if emperical cdf is 0 or 1;
data grad.logcurr5eth1;
	set grad.logcurr5eth;
	if ecdf_btc=1 or ecdf_eth=1 or ecdf_btc=0 or ecdf_eth=0 then delete;
run;

*Estimate dynamic t copula with AR(1);
proc mcmc data=grad.logcurr5eth1 nmc=2000 outpost=sim;
	parms lambda0 0.5 nu 0.3 sigma 0.2 alpha 0.5 beta 0.5;
	prior lambda0 nu alpha ~normal(0.4,sd=0.2);
	prior beta sigma~normal(0.5,sd=0.2);
	random lambda~normal(alpha+beta*lambda.l1,sd=sigma) subject=date icond=(lambda0) monitor=(lambda);
	phi=(1-exp(-lambda))/(1+exp(-lambda));
	ll=-1/2*log(abs(phi))+lgamma(1+2/nu)+lgamma(nu/2)+(1+nu/2)*log(1+1/(phi*nu)*(ecdf_btc**2+ecdf_eth**2))
		-lgamma((1+2/nu)**2)-(1/2+nu/2)*log(1+ecdf_btc**2/nu)-(1/2+nu/2)*log(1+ecdf_eth**2/nu);
	model ecdf_btc ecdf_eth ~ general(ll);
run;

