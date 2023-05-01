
// model in final form to be used on real data after tests on simulated data

functions {
  // define piecewise linear functions to describe the relationship between the 
  // parameters 
  real piecelm(real t, real tp, real wp, real wc, real x) {
    // dct linearly increases from 0 to x from the start to the peak of the infection
    if (t > (tp - wp) && t <= tp)
      return((x / wp) * (t + wp - tp));
    // dct linearly decreases from x to 0 from the peak to the end of the infection
    if (t > tp && t <= (tp + wc))
      return(x - ((x / wc) * (t - tp)));
    // dct is 0 if the time is outside of the infection episode 
    if (t <= (tp - wp) || t > (tp + wc) ) {
      return(0);
    }
    else 
      return(0);
  }
  
  real logdistsd(real varmean, real varsd) {
    return sqrt(log(1 + ((varsd / varmean) * (varsd / varmean))));
  }
  
  real logdistmean(real varmean, real logsd) {
    return log(varmean) - (0.5 * (logsd * logsd));
  }
  
}

data {
  int<lower = 0> N; // number of data points
  int<lower = 0> n_id; // number of infections (number of individuals for now, 
  // will need to represent number of infections in the future once including reinfections)
  real<lower = 0> lod; // limit of detection 
  int<lower = 0> id[N]; // vector for which individual the observation belongs too 
  real t[N]; // vector for the date of the observation
  real<lower = 0, upper = lod> ct[N]; // Ct value for a given observation
  real<lower = 0> dctpmean_p; // prior mean for peak dct (lod - ct)
  real<lower = 0> dctpsd_p; // prior sd for peak dct
  real<lower = 0> tpsd_p; // prior sd for date of peak of infection 
  real<lower = 0> wpmax; // prior max proliferation duration (onset to peak)
  real<lower = 0> wpmean_p; // prior mean proliferation duration
  real<lower = 0> wpsd_p; // prior sd proliferation duration 
  real<lower = 0> wcmax; // prior max clearance duration 
  real<lower = 0> wcmean_p; // prior mean clearance duration
  real<lower = 0> wcsd_p; // prior sd clearance duration
  real<lower = 0> sigma; // constant for sigma
}

transformed data {
  real<lower =0, upper = lod> dct[N]; // difference in lod and ct 

  for(i in 1:N) {
    dct[i] = lod - ct[i];
  }

}

parameters {
  real<lower=0, upper=lod> dctpmean; // population mean peak dct 
  real<lower=0, upper=wpmax> wpmean; // population mean proliferation duration
  real<lower=0, upper=wcmax> wcmean; // population mean clearance duration 

  real<lower=0> dctpsd; // population sd peak dct
  real<lower=0> wpsd; // population sd proliferation duration
  real<lower=0> wcsd; // population sd clearance duration

  real tp[n_id]; // peak date for each indivdiual
  real<lower=0, upper=lod> x[n_id];    // peak dct 
  real<lower=0, upper=wpmax> wp[n_id];  // proliferation duration
  real<lower=0, upper=wcmax> wc[n_id];  // clearance duration
 
}

transformed parameters {
  real mu[N];
  
  for(i in 1:N){
    mu[i]=piecelm(t[i], tp[id[i]], wp[id[i]], wc[id[i]], x[id[i]]);
  };
  
  real logdctpmean; 
  real logdctpsd;
  real logwpmean;
  real logwpsd;
  real logwcmean;
  real logwcsd;
  
  logdctpsd = logdistsd(dctpmean, dctpsd);
  logdctpmean = logdistmean(dctpmean, logdctpsd);
  logwpsd = logdistsd(wpmean, wpsd);
  logwpmean = logdistmean(wpmean, logwpsd);
  logwcsd = logdistsd(wcmean, wcsd);
  logwcmean = logdistmean(wcmean, logwcsd);
}

model {
  
  // Population mean/sd estimates
  
  dctpmean ~ normal(dctpmean_p, dctpsd_p) T[0, lod];
  wpmean ~ normal(wpmean_p, wpsd_p) T[0, wpmax];
  wcmean ~ normal(wcmean_p, wcsd_p) T[0, wcmax];

  // normal distributions for sd population parameters
  dctpsd ~ normal(dctpsd_p, dctpsd_p/2) T[0,];
  wpsd ~ normal(wpsd_p, wpsd_p/2) T[0,];
  wcsd ~ normal(wcsd_p, wcsd_p/2) T[0,];

  // individual parameters with log normal distribution - must scale the mean/sd parameters
  
  for(i in 1:n_id){
    tp[i] ~ normal(0,tpsd_p); 
    x[i] ~ lognormal(logdctpmean, logdctpsd) T[0, lod];
    wp[i] ~ lognormal(logwpmean, logwpsd) T[0, wpmax];
    wc[i] ~ lognormal(logwcmean, logwcsd) T[0, wcmax];
  }
  
  // Main model specification: 
  for(i in 1:N){

    dct[i] ~ normal(mu[i], sigma);
    if (dct[i] < 0 || dct[i] > lod)
      target += negative_infinity();
  
    else
      target += -log_diff_exp(
        normal_lcdf(lod | mu[i], sigma),
        normal_lcdf(0 | mu[i], sigma));
  
  }
  
}




