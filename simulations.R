library(doParallel)
cluster = makeCluster(detectCores())
clusterCall(cluster, function(x) .libPaths(x), .libPaths())
clusterCall(cluster, function(x) source("LDML.R"))
registerDoParallel(cluster)

source("LDML.R")

cate_fn <- function(x){(x[,1] + x[,2]) <= 0}
makedata.sim <- function(n,d=3,ovlp=3.5){
  x <- matrix(runif(n*d),n,d)*2-1
  p <- pnorm(ovlp*((x[,1])+(x[,3]))/2)
  t <- runif(n) <= p
  cate.true <- cate_fn(x)
  f0.true <- 0.
  f1.true <- f0.true + cate.true
  eps <- (1+x[,3])*rnorm(n)
  y0 <- f0.true + eps
  y1 <- f1.true + eps
  y <- (!t)*y0 + t*y1
  return(data.frame(x=x,t,y,p,y0,y1))
}

gammas = c(2/3)

d = 20
form_x = paste(paste("x.", 1:d, sep=""),collapse="+")
form_y = "y"
form_t = "t"
makedata =  function(n){makedata.sim(n, d, ovlp=3.0)}
ntest=1000000
datatest = makedata(ntest)
q.true = as.data.frame(do.call(rbind,lapply(gammas, function(gamma){data.frame(gamma=gamma,q1.true=quantile(datatest$y1, gamma),q0.true=quantile(datatest$y0, gamma))})))
d1.true = density(datatest$y1, n = 1, from = q.true$q1.true, to = q.true$q1.true, bw = 'SJ')$y
d0.true = density(datatest$y0, n = 1, from = q.true$q0.true, to = q.true$q0.true, bw = 'SJ')$y
rm(datatest)

K=5
trim=c(0.00001,0.99999)
trim.type='clip'
normalize=T
methods = list()
for(cls.meth.k in names(methods.classification)) {
  cls.meth = methods.classification[[cls.meth.k]]$method
  cls.opt  = methods.classification[[cls.meth.k]]$option
  methods[[paste('ipw' ,cls.meth.k)]] = list(method = est.quantile.ipw , option = list(method_prop=cls.meth, option_prop=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize, oracle.density=c(d0.true,d1.true)))
  methods[[paste('ldml',cls.meth.k)]] = list(method = est.quantile.ldml, option = list(method_ipw=cls.meth, option_ipw=cls.opt, method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize, oracle.density=c(d0.true,d1.true)))
  methods[[paste('dmlc',cls.meth.k)]] = list(method = est.quantile.dml , option = list(method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, cdf_regress=T, K=K, trim=trim, trim.type=trim.type, normalize=normalize, qrange=.01, oracle.density=c(d0.true,d1.true)))
  methods[[paste('dmlf',cls.meth.k)]] = list(method = est.quantile.dml , option = list(method_prop=cls.meth, option_prop=cls.opt, method_cdf=forestcdf, option_cdf=forestcdf_option, cdf_regress=F, K=K, trim=trim, trim.type=trim.type, normalize=normalize, qrange=.01, oracle.density=c(d0.true,d1.true)))
}
methods[['reg']] = list(method = est.quantile.dml , option = list(method_prop=NULL, option_prop=NULL, method_cdf=cls.meth, option_cdf=cls.opt, cdf_regress=T, K=K, trim=trim, trim.type=trim.type, normalize=normalize, qrange=.01, oracle.density=c(d0.true,d1.true)))


ns = c(100,200,400,800,1600,3200,6400,6400*2,6400*4)
seeds.data.per.run = 25
runs = 10
seeds.method = 1
methods_to_use = c(as.vector(outer(c('ipw','ldml','dmlc','dmlf'),c('forest'),paste)),'reg')

print(paste('processes per run',length(ns)*seeds.data.per.run*seeds.method*length(methods_to_use)))

results = list()

for (j in 1:runs) {
  print(paste('run',j,Sys.time()))
  res = foreach(n = ns, .combine='rbind', .inorder=FALSE)%:%foreach(seed.data = (1+(j-1)*seeds.data.per.run):(j*seeds.data.per.run), .combine='rbind', .inorder=FALSE)%:%foreach(seed.method = 1:seeds.method, .combine='rbind', .inorder=FALSE)%:%foreach(method = methods_to_use, .combine='rbind', .inorder=FALSE, .packages=packages.needed) %dopar% {
    set.seed(seed.data);
    data = makedata(n);
    set.seed(seed.method);
    ret=tryCatch(
      do.call(methods[[method]]$method, append(list(gammas=gammas, data=data, form_x=form_x, form_t=form_t, form_y=form_y), methods[[method]]$option))
      , error = function(err) data.frame(gamma=gammas, q1=NA, q0=NA, qte=NA, se1=NA, se0=NA, seqte=NA)) %>% mutate(
        n=n,
        seed.data=seed.data,
        seed.method=seed.method,
        method=method);
    rm(data);
    ret
  };
  results[[j]] = res
}

stopCluster(cluster)
results = do.call(rbind,results)

alpha = .1
zconf = qnorm(1-alpha/2)
g=gammas[[1]]
plotres = results %>%
  group_by(n, seed.data, gamma, method) %>% summarise(est=median(q1, na.rm = TRUE), se=median(se1, na.rm = TRUE), se.oracle=median(se1.oracle, na.rm = TRUE)) %>% 
  left_join(q.true, by = 'gamma') %>% mutate(sqerr = (est-q1.true)**2, err = est-q1.true, covered = abs(est-q1.true)<=zconf*se, covered.oracle = abs(est-q1.true)<=zconf*se.oracle)
plotres$method[plotres$method=='dmlc forest'] = 'DML-D'
plotres$method[plotres$method=='dmlf forest'] = 'DML-F'
plotres$method[plotres$method=='ipw forest'] = 'IPW'
plotres$method[plotres$method=='ldml forest'] = 'LDML'
plotres$method[plotres$method=='reg'] = 'PI'
plotres$method = factor(plotres$method, c('LDML','IPW','DML-D','DML-F','PI'))

library(psych)
plot_mse = plotres %>% ggplot() + aes(n, sqerr, linetype=method, shape=method, color=method, fill=method) + 
  stat_summary(fun = winsor.mean, geom = "line") + stat_summary(fun = winsor.mean, geom = "point") + 
  stat_summary(fun.data = ~list(ymin=winsor.mean(.)-winsor.sd(.)/sqrt(length(.)), ymax=winsor.mean(.)+winsor.sd(.)/sqrt(length(.))), geom = "ribbon", alpha = 0.1) + 
  scale_y_log10() + scale_x_log10()  + ylab('Mean-Squared Error') + xlab('n')
plot_coverage = plotres %>% mutate(covered=case_when(method!="PI"&method!="IPW"~covered,TRUE~NaN)) %>% group_by(method,n) %>% summarise(coverage = mean(covered), coveragese = sqrt(coverage*(1-coverage)/n())) %>%
  ggplot() + aes(n, coverage, ymax=coverage+coveragese, ymin=coverage-coveragese, linetype=method, shape=method, color=method, fill=method) + geom_line() + scale_x_log10() + geom_point() + geom_ribbon(alpha=.1) + ylab('Coverage') + xlab('n') + geom_segment(aes(x=min(plotres$n),xend=max(plotres$n),y=1-alpha,yend=1-alpha),color='black',linetype=3)
plot_coverage.oracle = plotres %>% mutate(covered=case_when(method!="PI"&method!="IPW"~covered.oracle,TRUE~NaN)) %>% group_by(method,n) %>% summarise(coverage = mean(covered), coveragese = sqrt(coverage*(1-coverage)/n())) %>%
  ggplot() + aes(n, coverage, ymax=coverage+coveragese, ymin=coverage-coveragese, linetype=method, shape=method, color=method, fill=method) + geom_line() + scale_x_log10() + geom_point() + geom_ribbon(alpha=.1) + ylab('Coverage') + xlab('n') + geom_segment(aes(x=min(plotres$n),xend=max(plotres$n),y=1-alpha,yend=1-alpha),color='black',linetype=3)

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ggsave('sims_leg.pdf',plot=g_legend(plot_mse))

ggsave('sims_mse.pdf',plot=plot_mse      + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)
ggsave('sims_cov.pdf',plot=plot_coverage + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)
ggsave('sims_cov_oracle.pdf',plot=plot_coverage.oracle + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)

ggsave('sims_mse_narrow.pdf',plot=plot_mse      + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 3.5)
ggsave('sims_cov_narrow.pdf',plot=plot_coverage + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 2.25)
ggsave('sims_cov_oracle_narrow.pdf',plot=plot_coverage.oracle + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 2.25)

ggsave('sims_mse_wide.pdf',plot=plot_mse    + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 5, width = 7)
ggsave('sims_cov_wide.pdf',plot=plot_coverage + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 5, width = 7)
ggsave('sims_cov_oracle_wide.pdf',plot=plot_coverage.oracle + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 5, width = 7)
