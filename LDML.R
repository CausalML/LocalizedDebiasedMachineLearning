packages.needed = c('foreach','tidyverse','Hmisc','gbm','glmnetUtils','nnet','hdm','ks','randomForest','quantregForest');
lapply(packages.needed, library, character.only = TRUE);

# returns list of length n of integers in {1,..,K} indicating fold membership
make.cvgroup = function(n, K, right = TRUE) {
  split     = runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right = right)))
}

# conditions on fraction of 1s to 0s being similar across folds
make.cvgroup.balanced = function(data, K, form_t) {
  cvgroup = numeric(nrow(data))
  cvgroup[data[[form_t]]==1] = make.cvgroup(sum(data[[form_t]]==1), K, right = TRUE)
  cvgroup[data[[form_t]]==0] = make.cvgroup(sum(data[[form_t]]==0), K, right = FALSE)
  return(cvgroup)
}

# conditions on distribution of (t,w) being similar across folds
make.cvgroup.balanced2 = function(data, K, form_t, form_w) {
  cvgroup = numeric(nrow(data))
  for (t in 0:1) { for (w in 0:1) {
    cvgroup[data[[form_t]]==t&data[[form_w]]==w] = make.cvgroup(sum(data[[form_t]]==t&data[[form_w]]==w), K, right = (t==1))
  } }
  return(cvgroup)
}

cross.fit.propensities = function(data, cvgroup, form_x, form_t, method_prop, option_prop, trim=c(0.01,0.99), trim.type='none', normalize=T, trainmask=T) {
  K = max(cvgroup)
  prop      = numeric(nrow(data))
  for (k in 1:K) {
    prop[cvgroup==k] = method_prop(data, (cvgroup!=k) & trainmask, cvgroup==k, form_x, form_t, option_prop)
  }
  if(trim.type == 'drop') {
    keep = (prop>trim[1] & prop<trim[2])
    prop[!keep] = 0.5
  } else {
    keep = rep(T, nrow(data))
  }
  if(trim.type == 'clip') {
    prop[prop<trim[1]] = trim[1]
    prop[prop>trim[2]] = trim[2]
  }
  if(normalize) {
    prop[keep] = data[[form_t]][keep]*prop[keep]*mean(data[[form_t]][keep]/prop[keep]) + (1.-data[[form_t]][keep])*(1.-(1.-prop[keep])*mean((1.-data[[form_t]][keep])/(1.-prop[keep])))
  }
  return(list(prop=prop, keep=keep))
}

## Return i for i such that summing v[1:i] is closest to c
## (If v is increasing, as it is for ipw and ldml, really should do this by golden section search)
solve.cumsum = function(v,c) {
  return(which.min(abs(cumsum(v)-c)))
}

## Estimate the density of data in X at point x with weights w
density_ = function(X, w, x) {
  if(all(w>=0)) {
    density(X, n = 1, from = x, to = x, weights = w/sum(w), bw = 'SJ')$y
  } else {
    kde(X, eval.points = x, w = w/sum(w)*length(w), binned = F)$estimate
  }
}

check.data = function(data, form_x, form_t, form_y) {
  stopifnot(
    all(sort(unique(data[[form_t]]))==c(F,T))
  )
}

check.data2 = function(data, form_x, form_w, form_t, form_y) {
  stopifnot(
    all(sort(unique(data[[form_w]]))==c(F,T))
  )
  check.data(data, form_x, form_t, form_y)
}

## This uses the estimating equation 1/n sum_i T_i I[Y_i<=theta] / e(X_i) = gamma
## Where e(x)=P(T=1|X=x)
## We cross-fit e using K folds to estimate e(X_i) by ehat_i
## We then solve the equation as the gamma quantile of the data {Y_i:T_i=1} reweighted by W_i=(n_1/n)/ehat_i
## Symmetrically for control outcome
## If avg.eqn is T we solve the average equation where each element is using out-of-fold nuisnaces
## If avg.eqn is F we solve the equation in each fold and then average the estimates
## oracle.density is an option to pass a guess for the density for use in the standard error estimation
est.quantile.ipw = function(gammas, data, form_x, form_t, form_y, method_prop, option_prop, K=5, trim=c(0.01,0.99), trim.type='none', normalize=T, avg.eqn=T, oracle.density=NULL) {
  data = data%>%arrange(!! sym(form_y))
  cvgroup   = make.cvgroup.balanced(data, K, form_t)
  prop = cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  W1 = prop$keep*data[[form_t]]/prop$prop
  W0 = prop$keep*(1.-data[[form_t]])/(1.-prop$prop)
  return(foreach(gamma=gammas, .combine=rbind)%do% {
    q1 = if(avg.eqn) data[[form_y]][solve.cumsum(W1/sum(prop$keep),gamma)] else foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W1*(cvgroup==k)/sum(prop$keep&(cvgroup==k)),gamma)]}/K;
    q0 = if(avg.eqn) data[[form_y]][solve.cumsum(W0/sum(prop$keep),gamma)] else foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W0*(cvgroup==k)/sum(prop$keep&(cvgroup==k)),gamma)]}/K;
    psi1 = (W1[prop$keep] * (data[[form_y]][prop$keep] <= q1) - gamma) / density_(data[[form_y]][data[[form_t]]==1 & prop$keep], 1./prop$prop[data[[form_t]]==1 & prop$keep], q1);
    psi0 = (W0[prop$keep] * (data[[form_y]][prop$keep] <= q0) - gamma) / density_(data[[form_y]][data[[form_t]]==0 & prop$keep], 1./(1.-prop$prop[data[[form_t]]==0 & prop$keep]), q0);
    se1 = sd(psi1) / sqrt(sum(prop$keep));
    se0 = sd(psi0) / sqrt(sum(prop$keep));
    seqte = sd(psi1-psi0) / sqrt(sum(prop$keep));
    if(is.null(oracle.density)){data.frame(
    gamma=gamma,
    q1=q1,
    q0=q0,
    qte=q1-q0,
    se1=se1,
    se0=se0,
    seqte=seqte)}else{
    psi1.oracle = (W1[prop$keep] * (data[[form_y]][prop$keep] <= q1) - gamma) / oracle.density[2];
    psi0.oracle = (W0[prop$keep] * (data[[form_y]][prop$keep] <= q0) - gamma) / oracle.density[1];
    se1.oracle = sd(psi1.oracle) / sqrt(sum(prop$keep));
    se0.oracle = sd(psi0.oracle) / sqrt(sum(prop$keep));
    seqte.oracle = sd(psi1.oracle-psi0.oracle) / sqrt(sum(prop$keep));
    data.frame(
    gamma=gamma,
    q1=q1,
    q0=q0,
    qte=q1-q0,
    se1=se1,
    se0=se0,
    seqte=seqte,
    se1.oracle=se1.oracle,
    se0.oracle=se0.oracle,
    seqte.oracle=seqte.oracle)}
    })
}

## This uses the estimating equation 1/n sum_i T_i I[Y_i<=theta] / e(X_i) = gamma - 1/n sum_i f(theta,X_i)*(1-T_i/e(X_i))
## Where f(theta,x)=P(Y<=theta|X=x,T=1)
## We use LDML cross-fitting with K folds
## Namely, for each fold, we take half the remaining data and use it for fitting an initial guess for theta
## This initial guess is done using cross-fit IPW as above on this data subset alone
## On the other half of the remaining data we fit e(.) and f(intialguess,.) and use these on X_i in the fold to get ehat_i, fhat_i
## (if semiadaptive is set to TRUE then we use all out-of-fold data for both IPW and fitting nuisances)
## We finally compute c1 = 1/n sum_i fhat_i*(1-T_i/ehat_i)
## And then solve the equation as the gamma-c1 quantile of the data {Y_i:T_i=1} reweighted by (n_1/n)/ehat_i
## Symmetrically for control outcome
## If q.oracle is given then we use that given fixed value for localization; it shuold be a data.frame with columns gamma, q1.true, q0.true
## If avg.eqn is T we solve the average equation where each element is using out-of-fold nuisnaces
## If avg.eqn is F we solve the equation in each fold and then average the estimates
## oracle.density is an option to pass a guess for the density for use in the standard error estimation
est.quantile.ldml = function(gammas, data, form_x, form_t, form_y, method_ipw, option_ipw, method_prop, option_prop, method_cdf, option_cdf, K=5, K_ipw=NULL, semiadaptive=FALSE, trim=c(0.01,0.99), trim.type='none', normalize=T, q.oracle=NULL, avg.eqn=T, oracle.density=NULL) {
  data = data%>%arrange(!! sym(form_y))
  cvgroup   = make.cvgroup.balanced(data, K, form_t)
  prop = cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  W1 = prop$keep*data[[form_t]]/prop$prop
  W0 = prop$keep*(1.-data[[form_t]])/(1.-prop$prop)
  if(is.null(K_ipw)) {K_ipw = ceil((K-1)/2)}
  if(is.null(q.oracle)) {
    ipwquant = foreach(k = 1:K, .combine=rbind)%do% {
    est.quantile.ipw(gammas, data[if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==0,], form_x, form_t, form_y, method_ipw, option_ipw, K=K_ipw, trim=trim, trim.type = trim.type, normalize = normalize) %>% mutate(k=k)
    }}
  return(foreach(gamma=gammas, .combine=rbind)%do% {
    cdf1      = numeric(nrow(data));
    cdf0      = numeric(nrow(data));
    for (k in 1:K) {
      ## take out k from the list of folds, renumber so k+1->k, k+2->k+1, ..., and use only the even folds after renumbering for ipw
      ## and use the odd folds after renumbering for fitting nuisances, using the result from ipq for eta_1
      ## unless semiadaptive is set to TRUE
      q1.ipw = if(is.null(q.oracle)) ipwquant%>%filter(gamma==!!gamma,k==!!k)%>%select(q1) else q.oracle%>%filter(gamma==!!gamma)%>%select(q1.true)
      q0.ipw = if(is.null(q.oracle)) ipwquant%>%filter(gamma==!!gamma,k==!!k)%>%select(q0) else q.oracle%>%filter(gamma==!!gamma)%>%select(q0.true)
      form_cdf1          = paste('I(',form_y,'<=',as.numeric(q1.ipw),')')
      cdf1[cvgroup==k]   = method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_t]]==1, cvgroup==k, form_x, form_cdf1, option_cdf)
      form_cdf0          = paste('I(',form_y,'<=',as.numeric(q0.ipw),')')
      cdf0[cvgroup==k]   = method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_t]]==0, cvgroup==k, form_x, form_cdf0, option_cdf)
    };
    q1 = if(avg.eqn) {
        data[[form_y]][solve.cumsum(W1/sum(prop$keep),gamma - mean(cdf1[prop$keep] * (1.- data[[form_t]][prop$keep]/prop$prop[prop$keep])))]
      } else {
        foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W1*(cvgroup==k)/sum(prop$keep&(cvgroup==k)),gamma - mean(cdf1[prop$keep&(cvgroup==k)] * (1.- data[[form_t]][prop$keep&(cvgroup==k)]/prop$prop[prop$keep&(cvgroup==k)])))]}/K};
    q0 = if(avg.eqn) {
        data[[form_y]][solve.cumsum(W0/sum(prop$keep),gamma - mean(cdf0[prop$keep] * (1.- (1.-data[[form_t]][prop$keep])/(1.-prop$prop[prop$keep]))))]
      } else {
        foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W0*(cvgroup==k)/sum(prop$keep&(cvgroup==k)),gamma - mean(cdf0[prop$keep&(cvgroup==k)] * (1.- (1.-data[[form_t]][prop$keep&(cvgroup==k)])/(1.-prop$prop[prop$keep&(cvgroup==k)]))))]}/K};
    psi1 = (W1[prop$keep] * (data[[form_y]][prop$keep] <= q1) - gamma - cdf1[prop$keep] * (1.- data[[form_t]][prop$keep]/prop$prop[prop$keep])) / density_(data[[form_y]][data[[form_t]]==1 & prop$keep], 1./prop$prop[data[[form_t]]==1 & prop$keep], q1);
    psi0 = (W0[prop$keep] * (data[[form_y]][prop$keep] <= q0) - gamma - cdf0[prop$keep] * (1.- (1.-data[[form_t]][prop$keep])/(1.-prop$prop[prop$keep]))) / density_(data[[form_y]][data[[form_t]]==0 & prop$keep], 1./(1.-prop$prop[data[[form_t]]==0 & prop$keep]), q0);
    se1 = sd(psi1) / sqrt(sum(prop$keep));
    se0 = sd(psi0) / sqrt(sum(prop$keep));
    seqte = sd(psi1-psi0) / sqrt(sum(prop$keep));
    if(is.null(oracle.density)){data.frame(
      gamma=gamma,
      q1 = q1,
      q0 = q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte
    )}else{
    psi1.oracle = (W1[prop$keep] * (data[[form_y]][prop$keep] <= q1) - gamma - cdf1[prop$keep] * (1.- data[[form_t]][prop$keep]/prop$prop[prop$keep])) / oracle.density[2];
    psi0.oracle = (W0[prop$keep] * (data[[form_y]][prop$keep] <= q0) - gamma - cdf0[prop$keep] * (1.- (1.-data[[form_t]][prop$keep])/(1.-prop$prop[prop$keep]))) / oracle.density[1];
    se1.oracle = sd(psi1.oracle) / sqrt(sum(prop$keep));
    se0.oracle = sd(psi0.oracle) / sqrt(sum(prop$keep));
    seqte.oracle = sd(psi1.oracle-psi0.oracle) / sqrt(sum(prop$keep));
    data.frame(
      gamma=gamma,
      q1 = q1,
      q0 = q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte,
      se1.oracle=se1.oracle,
      se0.oracle=se0.oracle,
      seqte.oracle=seqte.oracle
    )
    }
  })
}

## This uses the estimating equation 1/n sum_i T_i I[Y_i<=theta] / e(X_i) = gamma - 1/n sum_i f(theta,X_i)*(1-T_i/e(X_i))
## Where f(theta,x)=P(Y<=theta|X=x,T=1)
## We use non-localized DML cross-fitting with K folds
## Namely, for each fold, we take the remaining data and use it for fitting e(.) and the whole f(.,.)
## We use this on X_i in the fold to get ehat_i, fhat_i(theta)
## We restrict the range of theta to a discretized grid of marginal Y quantiles
## Then we fit f(theta,.) for each theta in the range
## The range of quantiles is given by the list qrange
## Or if qrange is a number the we use all quantiles by qrange increments
## We then solve the equation by brute force search over theta in qrange to minimize abs of eqn
## If cdf_regress is T then method_cdf is a binary regresison method that we apply to each I[Y_i<=theta]
## If cdf_regress is F then method_cdf takes list of quantiles to simultaneously predict
## If avg.eqn is T we solve the average equation where each element is using out-of-fold nuisnaces
## If avg.eqn is F we solve the equation in each fold and then average the estimates
## If method_prop is NULL then we remove the propensity correction, ie just invert the fitted CDF
## oracle.density is an option to pass a guess for the density for use in the standard error estimation
est.quantile.dml = function(gammas, data, form_x, form_t, form_y, method_prop, option_prop, method_cdf, option_cdf, cdf_regress=T, qrange=0.01, K=5, trim=c(0.01,0.99), trim.type='none', normalize=T, avg.eqn=T, oracle.density=NULL) {
  if(length(qrange)==1) {
    qrange = seq(qrange,1.-qrange,qrange)
    #qrange = qrange[(qrange>=min(gammas)-.1) & (qrange<=max(gammas)+.1)]
  }
  cvgroup   = make.cvgroup.balanced(data, K, form_t)
  if(is.null(method_prop)) {
    prop = list(prop=rep(0, nrow(data)), keep=rep(T, nrow(data)))
  } else {
    prop = cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  }
  yqs  = quantile(data[[form_y]], qrange)
  cdf1      = matrix(0L, length(qrange), nrow(data));
  cdf0      = matrix(0L, length(qrange), nrow(data));
  if (cdf_regress) {
    for (i in 1:length(qrange)) {
      for (k in 1:K) {
        form_cdf1          = paste('I(',form_y,'<=',as.numeric(yqs[i]),')')
        cdf1[i,cvgroup==k] = method_cdf(data, cvgroup!=k & data[[form_t]]==1, cvgroup==k, form_x, form_cdf1, option_cdf)
        form_cdf0          = paste('I(',form_y,'<=',as.numeric(yqs[i]),')')
        cdf0[i,cvgroup==k] = method_cdf(data, cvgroup!=k & data[[form_t]]==0, cvgroup==k, form_x, form_cdf0, option_cdf)
      }
    }
  } else {
    for (k in 1:K) {
      cdf1[,cvgroup==k] = t(method_cdf(data, cvgroup!=k & data[[form_t]]==1, cvgroup==k, form_x, form_y, yqs, option_cdf))
      cdf0[,cvgroup==k] = t(method_cdf(data, cvgroup!=k & data[[form_t]]==0, cvgroup==k, form_x, form_y, yqs, option_cdf))
    }
  }
  yleq = outer(yqs,data[[form_y]],'>=')
  if(is.null(method_prop)) {
    a1 = if(avg.eqn) (cdf1 %*% (prop$keep)) else foreach(k=1:K)%do%{cdf1 %*% (prop$keep&cvgroup==k)}
    a0 = if(avg.eqn) (cdf0 %*% (prop$keep)) else foreach(k=1:K)%do%{cdf0 %*% (prop$keep&cvgroup==k)}
  } else {
    a1 = if(avg.eqn) (yleq %*% (prop$keep*data[[form_t]]/prop$prop) + cdf1 %*% (prop$keep*(1-data[[form_t]]/prop$prop))) else foreach(k=1:K)%do%{(yleq %*% ((prop$keep&cvgroup==k)*data[[form_t]]/prop$prop) + cdf1 %*% ((prop$keep&cvgroup==k)*(1-data[[form_t]]/prop$prop)))}
    a0 = if(avg.eqn) (yleq %*% (prop$keep*(1-data[[form_t]])/(1-prop$prop)) + cdf0 %*% (prop$keep*(1-(1-data[[form_t]])/(1-prop$prop)))) else foreach(k=1:K)%do%{(yleq %*% ((prop$keep&cvgroup==k)*(1-data[[form_t]])/(1-prop$prop)) + cdf0 %*% ((prop$keep&cvgroup==k)*(1-(1-data[[form_t]])/(1-prop$prop))))}
  }
  return(foreach(gamma=gammas, .combine=rbind)%do% {
    q1 = if(avg.eqn) yqs[which.min(abs(a1/sum(prop$keep) - gamma))] else foreach(k=1:K, .combine=sum)%do%{yqs[which.min(abs(a1[[k]]/sum(prop$keep&cvgroup==k) - gamma))]}/K;
    q0 = if(avg.eqn) yqs[which.min(abs(a0/sum(prop$keep) - gamma))] else foreach(k=1:K, .combine=sum)%do%{yqs[which.min(abs(a0[[k]]/sum(prop$keep&cvgroup==k) - gamma))]}/K;
    i1 = which.min(abs(yqs-q1))
    i0 = which.min(abs(yqs-q0))
    if(is.null(method_prop)) {
      se1 = 0
      se0 = 0
      seqte = 0
    } else {
      psi1 = (yleq[i1,] * (prop$keep*data[[form_t]]/prop$prop) + cdf1[i1,] * (prop$keep*(1-data[[form_t]]/prop$prop)) - gamma) / density_(data[[form_y]][data[[form_t]]==1 & prop$keep], 1./prop$prop[data[[form_t]]==1 & prop$keep], q1);
      psi0 = (yleq[i0,] * (prop$keep*(1-data[[form_t]])/(1-prop$prop)) + cdf0[i0,] * (prop$keep*(1-(1-data[[form_t]])/(1-prop$prop))) - gamma) / density_(data[[form_y]][data[[form_t]]==0 & prop$keep], 1./(1.-prop$prop[data[[form_t]]==0 & prop$keep]), q0);
      se1 = sd(psi1[prop$keep]) / sqrt(sum(prop$keep));
      se0 = sd(psi0[prop$keep]) / sqrt(sum(prop$keep));
      seqte = sd(psi1[prop$keep]-psi0[prop$keep]) / sqrt(sum(prop$keep));
    }
    if(is.null(oracle.density)){data.frame(
      gamma=gamma,
      q1=q1,
      q0=q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte
    )}else{
    psi1.oracle = (yleq[i1,] * (prop$keep*data[[form_t]]/prop$prop) + cdf1[i1,] * (prop$keep*(1-data[[form_t]]/prop$prop)) - gamma) / oracle.density[2];
    psi0.oracle = (yleq[i0,] * (prop$keep*(1-data[[form_t]])/(1-prop$prop)) + cdf0[i0,] * (prop$keep*(1-(1-data[[form_t]])/(1-prop$prop))) - gamma) / oracle.density[1];
    se1.oracle = sd(psi1.oracle[prop$keep]) / sqrt(sum(prop$keep));
    se0.oracle = sd(psi0.oracle[prop$keep]) / sqrt(sum(prop$keep));
    seqte.oracle = sd(psi1.oracle[prop$keep]-psi0.oracle[prop$keep]) / sqrt(sum(prop$keep));
    data.frame(
      gamma=gamma,
      q1=q1,
      q0=q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte,
      se1.oracle=se1.oracle,
      se0.oracle=se0.oracle,
      seqte.oracle=seqte.oracle
    )}
  })
}

est.ivquantile.ipw = function(gammas, data, form_x, form_w, form_t, form_y, method_prop, option_prop, K=5, trim=c(0.01,0.99), trim.type='none', normalize=T, avg.eqn=T, one.way.noncompliance=F) {
  data = data%>%arrange(!! sym(form_y))
  cvgroup   = make.cvgroup.balanced2(data, K, form_t, form_w)
  propw   = cross.fit.propensities(data, cvgroup, form_x, form_w, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  proptw1 = cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trainmask = data[[form_w]]==1, trim=NULL, trim.type='none', normalize=F)
  proptw0 = if (one.way.noncompliance) NULL else cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trainmask = data[[form_w]]==0, trim=NULL, trim.type='none', normalize=F)
  keep    = propw$keep
  propw   = propw$prop
  proptw1 = proptw1$prop
  proptw0 = if (one.way.noncompliance) rep(0.,nrow(data)) else proptw0$prop
  nu1 = sum(keep*(proptw1 + data[[form_w]]*(data[[form_t]]-proptw1)/propw))/sum(keep)
  nu0 = if (one.way.noncompliance) 0. else sum(keep*(proptw0 + (1-data[[form_w]])*(data[[form_t]]-proptw0)/(1-propw)))/sum(keep)
  W  = keep*(data[[form_w]]-propw)/(propw*(1-propw))
  W1 = data[[form_t]]*W / (nu1 - nu0)
  W0 = (1-data[[form_t]])*W / (nu0 - nu1)
  return(foreach(gamma=gammas, .combine=rbind)%do% {
    q1 = if(avg.eqn) data[[form_y]][solve.cumsum(W1/sum(keep),gamma)] else foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W1*(cvgroup==k)/sum(keep&(cvgroup==k)),gamma)]}/K;
    q0 = if(avg.eqn) data[[form_y]][solve.cumsum(W0/sum(keep),gamma)] else foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W0*(cvgroup==k)/sum(keep&(cvgroup==k)),gamma)]}/K;
    psi1 = (W1[keep] * (data[[form_y]][keep] <= q1) - gamma) / (density_(data[[form_y]][data[[form_t]]==1 & keep], W[data[[form_t]]==1 & keep] / (nu1 - nu0), q1));
    psi0 = (W0[keep] * (data[[form_y]][keep] <= q0) - gamma) / (density_(data[[form_y]][data[[form_t]]==0 & keep], W[data[[form_t]]==0 & keep] / (nu0 - nu1), q0));
    se1 = sd(psi1) / sqrt(sum(keep));
    se0 = sd(psi0) / sqrt(sum(keep));
    seqte = sd(psi1-psi0) / sqrt(sum(keep));
    data.frame(
      gamma=gamma,
      q1=q1,
      q0=q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte
    )})
}

est.ivquantile.ldml = function(gammas, data, form_x, form_w, form_t, form_y, method_ipw, option_ipw, method_prop, option_prop, method_cdf, option_cdf, K=5, K_ipw=NULL, semiadaptive=FALSE, trim=c(0.01,0.99), trim.type='none', normalize=T, one.way.noncompliance=F, q.oracle=NULL, avg.eqn=T) {
  data = data%>%arrange(!! sym(form_y))
  cvgroup = make.cvgroup.balanced2(data, K, form_t, form_w)
  propw   = cross.fit.propensities(data, cvgroup, form_x, form_w, method_prop, option_prop, trim=trim, trim.type=trim.type, normalize=normalize)
  proptw1 = cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trainmask = data[[form_w]]==1, trim=NULL, trim.type='none', normalize=F)
  proptw0 = if (one.way.noncompliance) NULL else cross.fit.propensities(data, cvgroup, form_x, form_t, method_prop, option_prop, trainmask = data[[form_w]]==0, trim=NULL, trim.type='none', normalize=F)
  keep    = propw$keep
  propw   = propw$prop
  proptw1 = proptw1$prop
  proptw0 = if (one.way.noncompliance) rep(0.,nrow(data)) else proptw0$prop
  nu1 = sum(keep*(proptw1 + data[[form_w]]*(data[[form_t]]-proptw1)/propw))/sum(keep)
  nu0 = if (one.way.noncompliance) 0. else sum(keep*(proptw0 + (1-data[[form_w]])*(data[[form_t]]-proptw0)/(1-propw)))/sum(keep)
  W  = keep*(data[[form_w]]-propw)/(propw*(1-propw))
  W1 = data[[form_t]]*W / (nu1 - nu0)
  W0 = (1-data[[form_t]])*W / (nu0 - nu1)
  if(is.null(K_ipw)) {K_ipw = ceil((K-1)/2)}
  if(is.null(q.oracle)) {
    ipwquant = foreach(k = 1:K, .combine=rbind)%do% {
      est.ivquantile.ipw(gammas, data[if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==0,], form_x, form_w, form_t, form_y, method_ipw, option_ipw, K=K_ipw, trim=trim, trim.type = trim.type, normalize = normalize, one.way.noncompliance = one.way.noncompliance) %>% mutate(k=k)
    }}
  return(foreach(gamma=gammas, .combine=rbind)%do% {
    cdf11     = numeric(nrow(data)); # P(Y<=th,T=1|W=1,X)
    cdf10     = numeric(nrow(data)); # P(Y<=th,T=0|W=1,X)
    cdf01     = numeric(nrow(data)); # P(Y<=th,T=1|W=0,X)
    cdf00     = numeric(nrow(data)); # P(Y<=th,T=0|W=0,X)
    for (k in 1:K) {
      ## take out k from the list of folds, renumber so k+1->k, k+2->k+1, ..., and use only the even folds after renumbering for ipw
      ## and use the odd folds after renumbering for fitting nuisances, using the result from ipq for eta_1
      ## unless semiadaptive is set to TRUE
      q1.ipw = if(is.null(q.oracle)) ipwquant%>%filter(gamma==!!gamma,k==!!k)%>%select(q1) else q.oracle%>%filter(gamma==!!gamma)%>%select(q1.true)
      q0.ipw = if(is.null(q.oracle)) ipwquant%>%filter(gamma==!!gamma,k==!!k)%>%select(q0) else q.oracle%>%filter(gamma==!!gamma)%>%select(q0.true)
      form_cdf1           = paste('I(',form_y,'<=',as.numeric(q1.ipw),'&',form_t,'==1)')
      cdf11[cvgroup==k]   = method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_w]]==1, cvgroup==k, form_x, form_cdf1, option_cdf)
      cdf01[cvgroup==k]   = if (one.way.noncompliance) 0. else method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_w]]==0, cvgroup==k, form_x, form_cdf1, option_cdf)
      form_cdf0           = paste('I(',form_y,'<=',as.numeric(q1.ipw),'&',form_t,'==0)')
      cdf10[cvgroup==k]   = method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_w]]==1, cvgroup==k, form_x, form_cdf0, option_cdf)
      cdf00[cvgroup==k]   = method_cdf(data, (if(semiadaptive) cvgroup!=k else cvgroup!=k & (cvgroup-(cvgroup>k)) %% 2==1) & data[[form_w]]==0, cvgroup==k, form_x, form_cdf0, option_cdf)
    };
    q1 = if(avg.eqn) {
      data[[form_y]][solve.cumsum(W1/sum(keep), gamma - sum(keep * ( cdf11 - cdf01 - data[[form_w]] * cdf11 / propw + (1-data[[form_w]]) * cdf01 / (1-propw) ))/sum(keep)/(nu1 - nu0) )]
    } else {
      foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W1*(cvgroup==k)/sum(keep&(cvgroup==k)),gamma - sum((keep&cvgroup==k) ( cdf11 - cdf01 - data[[form_w]] * cdf11 / propw + (1-data[[form_w]]) * cdf01 / (1-propw) ))/sum(keep&cvgroup==k)/(nu1 - nu0) )]}/K};
    q0 = if(avg.eqn) {
      data[[form_y]][solve.cumsum(W0/sum(keep), gamma - sum(keep * ( cdf10 - cdf00 - data[[form_w]] * cdf10 / propw + (1-data[[form_w]]) * cdf00 / (1-propw) ))/sum(keep)/(nu0 - nu1) )]
    } else {
      foreach(k=1:K, .combine=sum)%do%{data[[form_y]][solve.cumsum(W0*(cvgroup==k)/sum(keep&(cvgroup==k)),gamma - (1/(nu0 - nu1)) * sum((keep&cvgroup==k) ( cdf10 - cdf00 - data[[form_w]] * cdf10 / propw + (1-data[[form_w]]) * cdf00 / (1-propw) ))/sum(keep&cvgroup==k) )]}/K};
    psi1 = (W1[keep] * (data[[form_y]][keep] <= q1) - gamma + ( cdf11[keep] - cdf01[keep] - data[[form_w]][keep] * cdf11[keep] / propw[keep] + (1-data[[form_w]][keep]) * cdf01[keep] / (1-propw[keep]) ) / (nu1 - nu0)) / (density_(data[[form_y]][data[[form_t]]==1 & keep], W[data[[form_t]]==1 & keep] / (nu1 - nu0), q1));
    psi0 = (W0[keep] * (data[[form_y]][keep] <= q0) - gamma + ( cdf10[keep] - cdf00[keep] - data[[form_w]][keep] * cdf10[keep] / propw[keep] + (1-data[[form_w]][keep]) * cdf00[keep] / (1-propw[keep]) ) / (nu0 - nu1)) / (density_(data[[form_y]][data[[form_t]]==0 & keep], W[data[[form_t]]==0 & keep / (nu0 - nu1)], q0));
    se1 = sd(psi1) / sqrt(sum(keep));
    se0 = sd(psi0) / sqrt(sum(keep));
    seqte = sd(psi1-psi0) / sqrt(sum(keep));
    data.frame(
      gamma=gamma,
      q1 = q1,
      q0 = q0,
      qte=q1-q0,
      se1=se1,
      se0=se0,
      seqte=seqte
    )
  })
}

const_option = list()
const = function(data, trainmask, testmask, form_x, form_resp, option) {
  rep(mean(if(form_resp%in%colnames(data)) data[[form_resp]] else model.matrix(as.formula(paste('~',form_resp,'-1')), data=data[trainmask,])[,2]), sum(testmask))
}

boost_option = list(distribution = 'bernoulli', bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=5, verbose = FALSE)
boost = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  fit = do.call(gbm, append(list(formula=form, data=data[trainmask,]), option));
  best = if('cv.folds' %in% names(option) && option[['cv.folds']]>0) gbm.perf(fit,plot.it=FALSE,method="cv") else gbm.perf(fit,plot.it=FALSE,method="OOB");
  return(predict(fit, n.trees=best, newdata=data[testmask,],  type="response"))
}

forest_option = list(nodesize=1, ntree=1000, na.action=na.omit, replace=TRUE)
forest = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste("as.factor(",form_resp, ") ~", form_x));
  tryCatch({
    fit = do.call(randomForest, append(list(formula=form, data=data[trainmask,]), option))
    return(predict(fit, newdata=data[testmask,],  type="prob")[,2])
  }, error = function(err) const(data, trainmask, testmask, form_x, form_resp, const_option))
}

neuralnet_option = list(linout=FALSE, size=2,  maxit=1000, decay=0.02, MaxNWts=10000,  trace=FALSE)
neuralnet = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  fit = do.call(nnet, append(list(formula=form, data=data[trainmask,]), option))
  return(predict(fit, newdata=data[testmask,],  type="raw"));
}

lassor_option = list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
lassor_post_option = append(lassor_option, list(post=TRUE))
lassor_nopost_option = append(lassor_option, list(post=FALSE))
lassor = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~ (", form_x,")^2"));
  tryCatch({
    fit = do.call(rlassologit, append(list(formula=form, data=data[trainmask,]), option))
    return(predict(fit, newdata=data[testmask,],  type="response"))
  }, error = function(err) const(data, trainmask, testmask, form_x, form_resp, const_option))
}

logistic_option = list(family = "binomial")
logistic = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  tryCatch({
    fit = do.call(glm, append(list(formula=form, data=data[trainmask,]), option))
    return(predict(fit, newdata=data[testmask,],  type="response"))
  }, error = function(err) const(data, trainmask, testmask, form_x, form_resp, const_option))
}

reglm_option = list(family="binomial")
reglm_lasso_option = append(reglm_option, list(alpha=1))
reglm_ridge_option = append(reglm_option, list(alpha=0))
reglm_elast_option = append(reglm_option, list(alpha=0.5))
reglm = function(data, trainmask, testmask, form_x, form_resp, option) {
  form = as.formula(paste(form_resp, "~ (", form_x,")^2"));
  tryCatch({
    fit = do.call(cv.glmnet, append(list(formula=form, data=data[trainmask,]), option))
    return(predict(fit, newdata=data[testmask,],  type="response"))
  }, error = function(err) const(data, trainmask, testmask, form_x, form_resp, const_option))
}

methods.classification = list(
  boost = list(method=boost, option=boost_option),
  forest = list(method=forest, option=forest_option),
  neuralnet = list(method=neuralnet, option=neuralnet_option),
  lassor = list(method=lassor, option=lassor_nopost_option),
  lassorpost = list(method=lassor, option=lassor_post_option),
  lasso = list(method=reglm, option=reglm_lasso_option),
  ridge = list(method=reglm, option=reglm_ridge_option),
  elast = list(method=reglm, option=reglm_elast_option),
  logistic = list(method=logistic, option=logistic_option),
  const = list(method=const, option=const_option)
)

forestcdf_option = list()
forestcdf = function(data, trainmask, testmask, form_x, form_resp, ths, option) {
  form = as.formula(paste(form_resp, "~", form_x));
  lmfit = lm(form,  x = TRUE, y = TRUE, data = data[trainmask,]);
  fit = do.call(quantregForest, append(list(x=lmfit$x[ ,-1], y=lmfit$y), option));
  return(predict(fit, newdata=data[testmask,], what=function(z){colMeans(outer(z,ths,'<='))}));
}
