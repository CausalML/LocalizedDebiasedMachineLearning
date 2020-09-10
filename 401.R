library(doParallel)
cluster = makeCluster(detectCores())
clusterCall(cluster, function(x) .libPaths(x), .libPaths())
clusterCall(cluster, function(x) source("LDML.R"))
registerDoParallel(cluster)

source("LDML.R")

library(xtable)
library(psych)
library(foreign)

data      = read.dta("sipp1991.dta");
form_y    = "net_tfa";
form_w    = "e401";
form_t    = "p401";
form_x    = "age + inc + educ + fsize + marr + twoearn + db + pira + hown" 
form_xl   = "poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown"
one.way.noncompliance = T

gammas = c(.25,.5,.75)

Ks = c(5,15,25)
trim = c(0.01, 0.99)
trim.type = 'clip'
normalize = T
cls.meth.ks = c('lassor','boost','neuralnet')
seeds.method = 100

res = foreach(K = Ks, .combine='rbind', .inorder=FALSE)%:%foreach(cls.meth.k = cls.meth.ks, .combine='rbind', .inorder=FALSE)%:%foreach(seed.method = 1:seeds.method, .combine='rbind', .inorder=FALSE, .packages=packages.needed) %dopar% {
  set.seed(seed.method)
  cls.meth = methods.classification[[cls.meth.k]]$method
  cls.opt = methods.classification[[cls.meth.k]]$option
  tryCatch(
    est.quantile.ldml(gammas, data, if(grepl('lasso',cls.meth.k,fixed = TRUE)) form_xl else form_x, form_w, form_y, method_ipw=cls.meth, option_ipw=cls.opt, method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize)
    , error = function(err) data.frame(gamma=gammas, q1=NA, q0=NA, qte=NA, se1=NA, se0=NA, seqte=NA)) %>% mutate(
      K=K,
      trim.type=trim.type,
      seed.method=seed.method,
      method=cls.meth.k)
}

res_l = foreach(K = Ks, .combine='rbind', .inorder=FALSE)%:%foreach(cls.meth.k = cls.meth.ks, .combine='rbind', .inorder=FALSE)%:%foreach(seed.method = 1:seeds.method, .combine='rbind', .inorder=FALSE, .packages=packages.needed) %dopar% {
  set.seed(seed.method)
  cls.meth = methods.classification[[cls.meth.k]]$method
  cls.opt = methods.classification[[cls.meth.k]]$option
  tryCatch(
    est.ivquantile.ldml(gammas, data, if(grepl('lasso',cls.meth.k,fixed = TRUE)) form_xl else form_x, form_w, form_t, form_y, method_ipw=cls.meth, option_ipw=cls.opt, method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize, one.way.noncompliance=one.way.noncompliance)
    , error = function(err) data.frame(gamma=gammas, q1=NA, q0=NA, qte=NA, se1=NA, se0=NA, seqte=NA)) %>% mutate(
      K=K,
      trim.type=trim.type,
      seed.method=seed.method,
      method=cls.meth.k)
}

res2 = res %>% group_by(gamma,K,trim.type,method) %>% summarise(across(c('q1','q0','qte'), list(mean=~winsor.mean(.,na.rm = TRUE),se=~winsor.sd(.,na.rm = TRUE)/sqrt(n()))), across(c('se1','se0','seqte'), ~winsor.mean(.,na.rm = TRUE))) %>% rename(sete=seqte)
col_names1 = c('gamma', 'K', 'trim.type', 'method')
col_names2 = c('est', 'se_split', 'se_lead', 'type')
res3 = do.call(rbind, lapply(c('1','0','te'), function(x){res2 %>% select(c(col_names1,starts_with(paste('q',x,sep='')),starts_with(paste('se',x,sep='')))) %>% mutate(type=x) %>% setNames(c(col_names1,col_names2))})) %>% mutate(se = sqrt(se_lead^2 + se_split^2))
print(xtable(res3 %>% ungroup %>% filter(type=='te') %>% mutate(est = sprintf('%.1f (%.1f)',est,se)) %>% select(gamma,K,method,est) %>% pivot_wider(names_from=method, values_from=est) %>% select(gamma,K,lassor,neuralnet,boost) %>% mutate(K=as.factor(K)) ),include.rownames=FALSE)

res2 = res_l %>% group_by(gamma,K,trim.type,method) %>% summarise(across(c('q1','q0','qte'), list(mean=~winsor.mean(.,na.rm = TRUE),se=~winsor.sd(.,na.rm = TRUE)/sqrt(n()))), across(c('se1','se0','seqte'), ~winsor.mean(.,na.rm = TRUE))) %>% rename(sete=seqte)
col_names1 = c('gamma', 'K', 'trim.type', 'method')
col_names2 = c('est', 'se_split', 'se_lead', 'type')
res3 = do.call(rbind, lapply(c('1','0','te'), function(x){res2 %>% select(c(col_names1,starts_with(paste('q',x,sep='')),starts_with(paste('se',x,sep='')))) %>% mutate(type=x) %>% setNames(c(col_names1,col_names2))})) %>% mutate(se = sqrt(se_lead^2 + se_split^2))
print(xtable(res3 %>% ungroup %>% filter(type=='te') %>% mutate(est = sprintf('%.1f (%.1f)',est,se)) %>% select(gamma,K,method,est) %>% pivot_wider(names_from=method, values_from=est) %>% select(gamma,K,lassor,neuralnet,boost) %>% mutate(K=as.factor(K)) ),include.rownames=FALSE)

quantSE =  function(x, p) {
  quant <- quantile(x,p)
  R <- sqrt(p*(1-p)/length(x))
  f <- sapply(quant, function(a,b) density(b,from=a,to=a,n=1)$y,x)
  tibble(p=p,quantile=quant,SE=R/f)
}
quantSE(data[[form_y]][data[[form_w]]==1], gammas) %>% left_join(quantSE(data[[form_y]][data[[form_w]]==0], gammas), by='p') %>% mutate(te=quantile.x-quantile.y, sete=sqrt(SE.x^2+SE.y^2)) %>% mutate(est = sprintf('%.1f (%.1f)',te,sete)) %>% select(p,est)
quantSE(data[[form_y]][data[[form_t]]==1], gammas) %>% left_join(quantSE(data[[form_y]][data[[form_t]]==0], gammas), by='p') %>% mutate(te=quantile.x-quantile.y, sete=sqrt(SE.x^2+SE.y^2)) %>% mutate(est = sprintf('%.1f (%.1f)',te,sete)) %>% select(p,est)

gammas2 = seq(.1,.9,.01)
trim = c(0.01, 0.99)
normalize = T
cls.meth.k = 'lassor'
K = 15
trim.type = 'clip'
seeds.method = 100

res_range = foreach(seed.method = 1:seeds.method, .combine='rbind', .inorder=FALSE, .packages=packages.needed) %dopar% {
  set.seed(seed.method)
  cls.meth = methods.classification[[cls.meth.k]]$method
  cls.opt = methods.classification[[cls.meth.k]]$option
  tryCatch(
    est.quantile.ldml(gammas2, data, if(grepl('lasso',cls.meth.k,fixed = TRUE)) form_xl else form_x, form_w, form_y, method_ipw=cls.meth, option_ipw=cls.opt, method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize)
    , error = function(err) data.frame(gamma=gammas, q1=NA, q0=NA, qte=NA, se1=NA, se0=NA, seqte=NA)) %>% mutate(
      K=K,
      trim.type=trim.type,
      seed.method=seed.method,
      method=cls.meth.k)
}

res_l_range = foreach(seed.method = 1:seeds.method, .combine='rbind', .inorder=FALSE, .packages=packages.needed) %dopar% {
  set.seed(seed.method)
  cls.meth = methods.classification[[cls.meth.k]]$method
  cls.opt = methods.classification[[cls.meth.k]]$option
  tryCatch(
    est.ivquantile.ldml(gammas2, data, if(grepl('lasso',cls.meth.k,fixed = TRUE)) form_xl else form_x, form_w, form_t, form_y, method_ipw=cls.meth, option_ipw=cls.opt, method_prop=cls.meth, option_prop=cls.opt, method_cdf=cls.meth, option_cdf=cls.opt, K=K, trim=trim, trim.type=trim.type, normalize=normalize, one.way.noncompliance=one.way.noncompliance)
    , error = function(err) data.frame(gamma=gammas, q1=NA, q0=NA, qte=NA, se1=NA, se0=NA, seqte=NA)) %>% mutate(
      K=K,
      trim.type=trim.type,
      seed.method=seed.method,
      method=cls.meth.k)
}

col_names1 = c('gamma', 'K', 'trim.type', 'method')
col_names2 = c('est', 'se_split', 'se_lead', 'type')

res_range2 = res_range %>% group_by(gamma,K,trim.type,method) %>% summarise(across(c('q1','q0','qte'), list(mean=~winsor.mean(.,na.rm = TRUE),se=~winsor.sd(.,na.rm = TRUE)/sqrt(n()))), across(c('se1','se0','seqte'), ~winsor.mean(.,na.rm = TRUE))) %>% rename(sete=seqte)
res_range3 = do.call(rbind, lapply(c('1','0','TE'), function(x){res_range2 %>% select(c(col_names1,starts_with(paste('q',x,sep='')),starts_with(paste('se',x,sep='')))) %>% mutate(type=x) %>% setNames(c(col_names1,col_names2))})) %>% mutate(se = sqrt(se_lead^2 + se_split^2))

res_l_range2 = res_l_range %>% group_by(gamma,K,trim.type,method) %>% summarise(across(c('q1','q0','qte'), list(mean=~winsor.mean(.,na.rm = TRUE),se=~winsor.sd(.,na.rm = TRUE)/sqrt(n()))), across(c('se1','se0','seqte'), ~winsor.mean(.,na.rm = TRUE))) %>% rename(sete=seqte)
res_l_range3 = do.call(rbind, lapply(c('1','0','TE'), function(x){res_l_range2 %>% select(c(col_names1,starts_with(paste('q',x,sep='')),starts_with(paste('se',x,sep='')))) %>% mutate(type=x) %>% setNames(c(col_names1,col_names2))})) %>% mutate(se = sqrt(se_lead^2 + se_split^2))

raw_quantiles = data %>% select(!!! rlang::syms(form_w), !!! rlang::syms(form_y)) %>% group_by(!!! rlang::syms(form_w)) %>% summarise(q = quantile(!!! rlang::syms(form_y), gammas2), gamma=gammas2) %>% rename(type = (!!! rlang::syms(form_w)))
raw_quantiles = raw_quantiles %>% filter(type == 1) %>% left_join(raw_quantiles %>% filter(type == 0), by = 'gamma') %>% mutate(q = q.x-q.y, type='TE') %>% select(type,q,gamma) %>% bind_rows(.,raw_quantiles %>% mutate(type = as.factor(type)))

raw_l_quantiles = data %>% select(!!! rlang::syms(form_t), !!! rlang::syms(form_y)) %>% group_by(!!! rlang::syms(form_t)) %>% summarise(q = quantile(!!! rlang::syms(form_y), gammas2), gamma=gammas2) %>% rename(type = (!!! rlang::syms(form_t)))
raw_l_quantiles = raw_quantiles %>% filter(type == 1) %>% left_join(raw_quantiles %>% filter(type == 0), by = 'gamma') %>% mutate(q = q.x-q.y, type='TE') %>% select(type,q,gamma) %>% bind_rows(.,raw_quantiles %>% mutate(type = as.factor(type)))

combined   = rbind(raw_quantiles %>% rename(est=q,target=type) %>% mutate(se=NA,type='raw'), res_range3 %>% ungroup %>% select(gamma,est,se,type) %>% rename(target=type) %>% mutate(type='LDML'))
combined_l = rbind(raw_l_quantiles %>% rename(est=q,target=type) %>% mutate(se=NA,type='raw'), res_l_range3 %>% ungroup %>% select(gamma,est,se,type) %>% rename(target=type) %>% mutate(type='LDML'))

confidence = 0.9
zconf = qnorm(1-(1-confidence)/2)

plot_range   = combined   %>% mutate(target = case_when(target=='1'~'Eligible',target=='0'~'Ineligible',TRUE~'Effect'))        %>% mutate(target = factor(target,levels=c('Effect','Ineligible','Eligible'))) %>% ggplot() + aes(x=gamma,y=est,ymin=est-zconf*se,ymax=est+zconf*se, color=target, linetype=type, fill=target) + geom_line() + geom_ribbon(alpha = 0.25, color=NA) + xlab('Quantile') + ylab('Dollars')
plot_range_l = combined_l %>% mutate(target = case_when(target=='1'~'Participant',target=='0'~'Nonparticipant',TRUE~'Effect')) %>% ggplot() + aes(x=gamma,y=est,ymin=est-zconf*se,ymax=est+zconf*se, color=target, linetype=type, fill=target) + geom_line() + geom_ribbon(alpha = 0.25, color=NA) + xlab('Quantile') + ylab('Dollars')

ggsave('401kQTE.pdf',plot=plot_range    + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)
ggsave('401kLQTE.pdf',plot=plot_range_l + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)

stopCluster(cluster)
