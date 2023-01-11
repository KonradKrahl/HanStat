#' @title LinReg
#' @description A simple multiple linear regression function (OLS) and it's requirements
#' @param dv dependent vaiable name as a string
#' @param iv a string vector with the names of the independent variables, separated by commas, use c(iv_1,iv_2...iv_n)
#' @param data a data frame containing the variables
#' @param bootstrap set bootstap to TRUE or FALSE, if FALSE Number of bootstraps are ignored
#' @param outlier_controll set outlier_controll to TRUE or FALSE, to use cooks distance to exlude outliers, if bootstrap==TRUE, outlier_controll must be FALSE
#' @param set plot to TRUE to create simple scatterplots of correlation between variables
#' @return the results of linear regression, plots and all requirements plus an conclusion about the violations
#' @keywords Linear Regression
#' @export
#' @examples
#' LinReg(dependent_variable,c('independent_variable_1','independent_variable_2','independent_variable_3'),dataframe, bootstrap = TRUE, Number_Bootstrapps=1000, outlier_controll = FALSE, plot=TRUE)






LinReg<-function(dv,iv,data,bootstrap,Number_Bootstrapps, outlier_controll, plot){
  ##intstall if required
  list.of.packages <- c('car','lmtest','olsrr','boot','ggplot2','gridExtra')
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  ###packages
  suppressMessages(
    for (pkg in c('car','lmtest','olsrr','boot','ggplot2','gridExtra','crayon'))
    {
      eval(bquote(library(.(pkg))))
    }
  )
  ###data
  y<-data[dv]
  x<-data[iv]
  colnames(x)<-iv
  ###as numeric
  yx<-cbind(y,x)
  for (i in c(1:length(yx)))
    {
    yx[,i]<-as.numeric(yx[,i])
  }
  ###Cooks Distance and BS Conditions
  if (outlier_controll==TRUE & bootstrap==TRUE)
  {
    cat(red('->->->Chose outlier controll or Bootstrapping! Not both!'))
  }
  else if(outlier_controll==FALSE & bootstrap==FALSE)
  {
  ###assign dta frame object to environment
  .GlobalEnv$data_regression <- yx
  ###model
  ####no columns
  n<-ncol(data_regression)
  ####no of i variables
  Variables <- sprintf("data_regression[,%i]", c(1:n))
  string_list<-c()
  string_list[[1]]<-paste(Variables[1],'~')
  for (i in c(2:length(Variables)))
    {
    variable<-ifelse(i!=length(Variables), paste(Variables[i],'+'), Variables[i])
    string_list[[i]]<-variable
    }
  Variables<-unlist(string_list)
  formula_model<-paste(Variables,collapse = '')
  ####model
  print('-------------------------------------------------------------------------')
  print('---------------------------LINEAR MODEL RESULTS START--------------------')
  print('-------------------------------------------------------------------------')
  print('                                                                         ')
  print('                                                                         ')
  model<-lm(formula_model)
  names(model$coefficients)<-c('Intercept',iv)
  print(summary(model))
  print('---------------------------MODEL INTERPRETATION START--------------------')
  print('                                                                         ')
  print('                                                                         ')
  Fstat<-round(summary(model)$fstatistic[1],digits=3)
  df1<-round(summary(model)$fstatistic[2],digits=3)
  df2<-round(summary(model)$fstatistic[3],digits=3)
  p<-round(pf(Fstat,df1,df2,lower.tail=F),digits=3)
  r2<-round(summary(model)$adj.r.squared,digits=3)
  effect<-ifelse(r2/(1-r2)<.02,'no effect',
                 ifelse(r2/(1-r2)<.15 & r2/(1-r2)>=.02, print('small effect'),
                        ifelse(r2/(1-r2)<.35 & r2/(1-r2)>=.15, print('medium effect'),
                               'large effect')))
  ifelse(p<.05, print(paste('The conducted model can explain a significant proportion of variance of',dv,'by the included independent variables (F(',df1,',',df2,')=',Fstat,', R²=',r2,', p=',p,').')),
                print(paste('The conducted model can not explain a significant proportion of variance of',dv,'by the included independent variables (F(',df1,',',df2,')=',Fstat,', R²=',r2,', p=',p,').')))
  print(paste('According to Cohen(1992) this result can be described as a',effect,'(f²=',round(r2/(1-r2),digits=3),').'))
  print('                                                                         ')
  print('                                                                         ')
  print('---------------------------MODEL INTERPRETATION FIN----------------------')
  print('                                                                         ')
  print('                                                                         ')
  print('---------------------COEFFICIENTS INTERPRETATION START-------------------')
  print('                                                                         ')
  print('                                                                         ')
  for (i in c(1:length(iv)))
    {
    print(paste('An increase in',iv[i],'by one unit results in an increase in',dv,'by',round(model$coefficients[i+1],digits=3)))
    ifelse(summary(model)$coefficients[i,4]<.05, print('***The result is significant***'),print('The result is not significant'))
    print('                                                                       ')
    print('                                                                       ')
  }
  print('---------------------COEFFICIENTS INTERPRETATION FIN---------------------')
  print('-------------------------------------------------------------------------')
  print('------------------------LINEAR MODEL RESULTS FIN-------------------------')
  print('                                                                         ')
  print('                                                                         ')
  ###requiremnts ifelese
  print('-----------------------------REQUIREMENTS:-------------------------------')
  print('                                                                         ')
  print('                                                                         ')
  No_Violations=0
  ##DW
  print('----------------------------Autocorrelation------------------------------')
  print('                                                                         ')
  x<-durbinWatsonTest(model, max.lag=1, simulate=TRUE, reps=1000,
                     method=c("normal"),
                     alternative=c("two.sided"))
  print(paste('Durbin-Watson-Statistics: ',x[2]))
  print(paste('Autocorrelation: ',x[1]))
  if (x[3]>0.05)
    {
      print('No evidence for autocorrelation')
    } else {
      print('Autocorrelation detected!')
      No_Violations=No_Violations+1
    }
  print('                                                                         ')
  ##Linearity
  print('--------------------------------Liniarity--------------------------------')
  print('                                                                         ')
  x<-resettest(model, power = 2:3, type = c("fitted", "regressor",
                                               "princomp"), data = list())
  print(paste('RESET Statistic: ', x$statistic))
  if (x$p.value>0.05)
    {
      print('No evidence for non linearity')
    } else {
      print('evidence for non linearity!')
      No_Violations=No_Violations+1
    }
  print('                                                                         ')
  ##Multicoll
  print('----------------------------Multicollinearity----------------------------')
  print('                                                                         ')
  x<-ols_vif_tol(model)
  print(paste('Variance inflation factors: ', x$VIF))
  if (any(x$VIF>5)==FALSE)
    {
      print('No evidence for Multicollinearity')
    } else {
      print('Evidence for Multicollinearity!')
      No_Violations=No_Violations+1
    }
  print('                                                                         ')
  ##Homoscedasticity
  print('-----------------------------Homoscedasticity----------------------------')
  print('                                                                         ')
  x<-bptest(model, studentize=FALSE)
  print(paste('Breusch Pagan Test :', x$statistic))
  if (x$p.value>0.05)
    {
      print('No evidence for Heteroscedasticity')
    } else {
      print('Evidence for Heteroscedasticity!')
      No_Violations=No_Violations+1
    }
  print('                                                                         ')
  #Normal dependent
  print('--------------Residual normality distribution dependent variable---------')
  print('                                                                         ')
  Residuen<-resid(model)
  x<-shapiro.test(Residuen)
  print(paste('Shapiro-wilk Statistic :', x$statistic))
  if (x$p.value>0.05)
    {
      print('No evidence for violation of Gaussian assumption')
    } else {
      print('Evidence for violation of Gaussian assumption!')
      No_Violations=No_Violations+1
    }
  print('                                                                         ')
  print('                                                                         ')
  if (No_Violations==0)
  {
    cat(green(paste('->->->','Total number of violated requirements:',No_Violations,'!')))
  }else{
    cat(red(paste('->->->','Total number of violated requirements:',No_Violations,'!')))
  }
  ###outliers
  }
  else if (outlier_controll==TRUE & bootstrap==FALSE){
  ###assign dta frame object to environment
  .GlobalEnv$data_regression <- yx
  ###model
  ####no columns
  n<-ncol(data_regression)
  ####no of i variables
  Variables <- sprintf("data_regression[,%i]", c(1:n))
  string_list<-c()
  string_list[[1]]<-paste(Variables[1],'~')
  for (i in c(2:length(Variables)))
    {
    variable<-ifelse(i!=length(Variables), paste(Variables[i],'+'), Variables[i])
    string_list[[i]]<-variable
    }
  Variables<-unlist(string_list)
  formula_model<-paste(Variables,collapse = '')
  ###Cooks Distance
  model_raw<-lm(formula_model)
  CooksD<-cooks.distance(model_raw)
  yx$Valid<-(CooksD<(4/(length(yx[,1]))))
  yx<-subset(yx,yx$Valid==TRUE)
  yx$Valid<-NULL
  .GlobalEnv$data_regression <- yx
  ####model
  print('-------------------------------------------------------------------------')
  print('---------------------------LINEAR MODEL RESULTS START--------------------')
  print('-------------------------------------------------------------------------')
  print('                                                                         ')
  print('                                                                         ')
  model<-lm(formula_model)
  names(model$coefficients)<-c('Intercept',iv)
  print(summary(model))
  print('---------------------------MODEL INTERPRETATION START--------------------')
  print('                                                                         ')
  print('                                                                         ')
  Fstat<-round(summary(model)$fstatistic[1],digits=3)
  df1<-round(summary(model)$fstatistic[2],digits=3)
  df2<-round(summary(model)$fstatistic[3],digits=3)
  p<-round(pf(Fstat,df1,df2,lower.tail=F),digits=3)
  r2<-round(summary(model)$adj.r.squared,digits=3)
  effect<-ifelse(r2/(1-r2)<.02,'no effect',
                 ifelse(r2/(1-r2)<.15 & r2/(1-r2)>=.02, print('small effect'),
                        ifelse(r2/(1-r2)<.35 & r2/(1-r2)>=.15, print('medium effect'),
                               'large effect')))
  ifelse(p<.05, print(paste('The conducted model can explain a significant proportion of variance of',dv,'by the included independent variables (F(',df1,',',df2,')=',Fstat,', R²=',r2,', p=',p,').')),
         print(paste('The conducted model can not explain a significant proportion of variance of',dv,'by the included independent variables (F(',df1,',',df2,')=',Fstat,', R²=',r2,', p=',p,').')))
  print(paste('According to Cohen(1992) this result can be described as a',effect,'(f²=',round(r2/(1-r2),digits=3),').'))
  print('                                                                         ')
  print('                                                                         ')
  print('---------------------------MODEL INTERPRETATION FIN----------------------')
  print('                                                                         ')
  print('                                                                         ')
  print('---------------------COEFFICIENTS INTERPRETATION START-------------------')
  print('                                                                         ')
  print('                                                                         ')
  for (i in c(1:length(iv)))
  {
    print(paste('An increase in',iv[i],'by one unit results in an increase in',dv,'by',round(model$coefficients[i+1],digits=3)))
    ifelse(summary(model)$coefficients[i,4]<.05, print('***The result is significant***'),print('The result is not significant'))
    print('                                                                       ')
    print('                                                                       ')
  }
  print('---------------------COEFFICIENTS INTERPRETATION FIN---------------------')
  print('-------------------------------------------------------------------------')
  print('------------------------LINEAR MODEL RESULTS FIN-------------------------')
  print('                                                                         ')
  print('                                                                         ')
  ###requiremnts ifelese
  print('-----------------------------REQUIREMENTS:-------------------------------')
  print('                                                                         ')
  print('                                                                         ')
  No_Violations=0
  ##DW
  print('----------------------------Autocorrelation------------------------------')
  print('                                                                         ')
  x<-durbinWatsonTest(model, max.lag=1, simulate=TRUE, reps=1000,
                      method=c("normal"),
                      alternative=c("two.sided"))
  print(paste('Durbin-Watson-Statistics: ',x[2]))
  print(paste('Autocorrelation: ',x[1]))
  if (x[3]>0.05)
  {
    print('No evidence for autocorrelation')
  } else {
    print('Autocorrelation detected!')
    No_Violations=No_Violations+1
  }
  print('                                                                         ')
  ##Linearity
  print('--------------------------------Liniarity--------------------------------')
  print('                                                                         ')
  x<-resettest(model, power = 2:3, type = c("fitted", "regressor",
                                            "princomp"), data = list())
  print(paste('RESET Statistic: ', x$statistic))
  if (x$p.value>0.05)
  {
    print('No evidence for non linearity')
  } else {
    print('evidence for non linearity!')
    No_Violations=No_Violations+1
  }
  print('                                                                         ')
  ##Multicoll
  print('----------------------------Multicollinearity----------------------------')
  print('                                                                         ')
  x<-ols_vif_tol(model)
  print(paste('Variance inflation factors: ', x$VIF))
  if (any(x$VIF>5)==FALSE)
  {
    print('No evidence for Multicollinearity')
  } else {
    print('Evidence for Multicollinearity!')
    No_Violations=No_Violations+1
  }
  print('                                                                         ')
  ##Homoscedasticity
  print('-----------------------------Homoscedasticity----------------------------')
  print('                                                                         ')
  x<-bptest(model, studentize=FALSE)
  print(paste('Breusch Pagan Test :', x$statistic))
  if (x$p.value>0.05)
  {
    print('No evidence for Heteroscedasticity')
  } else {
    print('Evidence for Heteroscedasticity!')
    No_Violations=No_Violations+1
  }
  print('                                                                         ')
  #Normal dependent
  print('--------------Residual normality distribution dependent variable---------')
  print('                                                                         ')
  Residuen<-resid(model)
  x<-shapiro.test(Residuen)
  print(paste('Shapiro-wilk Statistic :', x$statistic))
  if (x$p.value>0.05)
  {
    print('No evidence for violation of Gaussian assumption')
  } else {
    print('Evidence for violation of Gaussian assumption!')
    No_Violations=No_Violations+1
  }
  print('                                                                         ')
  print('                                                                         ')
  if (No_Violations==0)
  {
    cat(green(paste('->->->','Total number of violated requirements:',No_Violations,'!')))
  }else{
    cat(red(paste('->->->','Total number of violated requirements:',No_Violations,'!')))
  }
  ###bootstrap
  }
  else{
  .GlobalEnv$data_regression <- yx
  ###formula
  ind_var=iv[1]
  if (length(iv)>1){
    for (i in c(2:length(iv))){
    ind_var=paste(ind_var,'+',iv[i])
    }
  }
  bs_formula=paste(dv,'~',ind_var)
  # Bootstrap 95% CI for R-Squared
  # function to obtain R-Squared from the data
    rsq <- function(formula, data, indices) {
    d <- data[indices,] # allows boot to select sample
    fit <- lm(formula, data=d)
    return(summary(fit)$r.square)
    }
  results <- boot(data=data_regression, statistic=rsq,
                    R=Number_Bootstrapps, formula=bs_formula)
  print(paste('R-Squared: ', results[1]))
  print(paste('R-Squared CI-95%: ', results[1]))
  print(results)
  print(boot.ci(results, type="bca"))
  # Bootstrap 95% CI for regression coefficients
  # function to obtain regression weights
  bs <- function(formula, data, indices) {
    d <- data[indices,] # allows boot to select sample
    fit <- lm(formula, data=d)
    return(coef(fit))
  }
  # bootstrapping with 1000 replications
  results <- boot(data=data_regression, statistic=bs,
                  R=Number_Bootstrapps, formula=bs_formula)
  # view results
  print(results)
  #####!!!For loop for no index)
  # get 95% confidence intervals
  print(boot.ci(results, type="bca", index=1)) # intercept)
  print(boot.ci(results, type="bca", index=2)) # wt
  print(boot.ci(results, type="bca", index=3)) # disp
  }
  ###plot
  statement_out_boot<-ifelse(bootstrap==TRUE & outlier_controll==TRUE,1,0)
  if (plot==TRUE & statement_out_boot==0)
  {
  ###+adjust grid.arrange
  plots<-c()
  counter=1
  for (i in c(2:(length(iv)+length(dv))))
  {
  yy<-colnames(data_regression)[1]
  xx<-colnames(data_regression)[i]
  plots[[counter]]<-eval(substitute(ggplot(data_regression, aes(x=data_regression[,(i)], y=data_regression[,1])) +
  geom_point()+
  geom_smooth(method='lm', linetype="dotted",color="darkred")+
  labs(title=paste('Relationship between ',xx,' & ',yy), x=xx, y=yy),list(i = i)))
  counter=counter+1
  }
  do.call("grid.arrange", c(plots, ncol = ceiling(length(iv)/2), nrow = ceiling(length(iv)/2)))
  }
}


