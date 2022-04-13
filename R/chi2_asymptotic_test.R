#' Chi2 based power analysis.
#' @description  Genotypes frequencies are computed separately for cases and
#' controls based on the prevalence of the diseases, the minor allele frequency of the 2 snps and their ORs.
#' n_cases cases and n_controls are then simulated using the genotype frequencies and used in a logistic regression
#' model to predict the diseases phenotype. Under the null hypothesis, the Wald test statistics of the coefficients follows a Chi-square distribution with 1df.
#' Under the alternative hypothesis, the Wald test statistic for the coefficients
#' follow a chi-square distribution with a non-centrality parameter equal to its expectation. The power is computed as the cumulative distribution of the
#' Wald test-static under the alternative hypothesis up to the 1-alpha quantile of the Wald-test statistic under the null.
#'
#' @param n_cases Number of cases
#' @param n_controls Number of controls
#' @param maf_snp1 maf of the snp1
#' @param maf_snp2 maf of the snp2
#' @param a_OR Odds Ratio for the additive effect of snp1
#' @param b_OR Odds Ratio for the additive effect of snp2
#' @param int_OR Odds Ratio for the interaction effect of snp1 and snp2
#' @param prevalence Prevalence of the disease
#' @param penetrance_ref_genotype Penetrance of the reference genotype
#' @param alpha the alpha threshold for overall significance
#' @param encoding the encoding of snp1 and snp2. Possibilities are "dosage","dominant","recessive","hetero","dosage_recessive","dominant_"recessive","hetero_recessive".
#' @param input_type the type of input; either "OR" for Odds ratio or "GRR" for Genotypic Risk Ratio
#'
#' @return a sinlge float. The power
#' @export
#'
#' @examples
#'
#' compute_chi2_power(n_cases = 1000,n_controls = 1000,maf_snp1 = 0.1,maf_snp2 = 0.1,
#' a_OR = 1,b_OR = 2,int_OR = 1.2,penetrance_ref_genotype = 0.1,
#' alpha = 1e-5,encoding = "dosage",input_type="OR")



compute_chi2_power=function(n_cases,n_controls,maf_snp1,maf_snp2,a_OR,b_OR,int_OR,a1,a2,b1,b2,i,prevalence,penetrance_ref_genotype,alpha,encoding="dosage",input_type){
  if (input_type=="OR"){
  stratified_geno_proba=geno_probas_from_OR(penetrance_ref_genotype=penetrance_ref_genotype,mafa=maf_snp1,mafb= maf_snp2,a_OR = a_OR,b_OR=b_OR,int_OR=int_OR,encoding=encoding)
  }
  else if (input_type=="GRR"){
    stratified_geno_proba=geno_probas_from_GRR(prevalence=prevalence,mafa=maf_snp1,mafb=maf_snp2,a1=a1,a2=a2,b1=b1,b2=b2,i=i)
  }
  cases33 = round(stratified_geno_proba$cases_geno_probs*n_cases,digits=0)
  controls33 = round(stratified_geno_proba$controls_geno_probs*n_controls,digits=0)
  d = as.numeric(t(cases33/(cases33+controls33))) # this is the proportion of cases with that genotype -> the outcome of the logistic regression
  weig = as.numeric(t(cases33+controls33))
  # the model encoding takes implicitly place here :
  if (encoding=="dosage"){
    aA = c( 0,0,0,
            1,1,1,
            2,2,2)
    aB = c(0,1,2,
            0,1,2,
           0,1,2)
  } else if(encoding=="dominant"){
    aA = c( 0,0,0,
            1,1,1,
            1,1,1)
    aB = c(0,1,1,
           0,1,1,
           0,1,1)

  }else if(encoding=="recessive"){
    aA = c( 0,0,0,
            0,0,0,
            1,1,1)
    aB = c(0,0,1,
           0,0,1,
           0,0,1)

  }else if(encoding=="hetero"){
    aA = c( 0,0,0,
            1,1,1,
            0,0,0)
    aB = c(0,1,0,
           0,1,0,
           0,1,0)

  }else if(encoding=="dosage_recessive"){
    aA = c( 0,0,0,
            1,1,1,
            2,2,2)
    aB = c(0,0,1,
           0,0,1,
           0,0,1)
  }else if(encoding=="dominant_recessive"){
    aA = c( 0,0,0,
            1,1,1,
            1,1,1)
    aB = c(0,0,1,
           0,0,1,
           0,0,1)
  }else if(encoding=="hetero_recessive"){
    aA = c( 0,0,0,
            1,1,1,
            0,0,0)
    aB = c(0,0,1,
           0,0,1,
           0,0,1)
  }
  epi.log = stats::glm( d ~ aA + aB + aA*aB, weights = weig , family = "binomial")
  if (nrow(summary(epi.log)$coef)==4)
  {
    int_est = exp(summary(epi.log)$coef["aA:aB","Estimate"])
    int_pval = summary(epi.log)$coef["aA:aB","Pr(>|z|)"]

    a_est=exp(summary(epi.log)$coef["aA","Estimate"])
    a_pval= summary(epi.log)$coef["aA","Pr(>|z|)"]

    b_est=exp(summary(epi.log)$coef["aB","Estimate"])
    b_pval= summary(epi.log)$coef["aB","Pr(>|z|)"]

    # Compute the power
    q.thres = stats::qchisq(alpha,df=1,lower.tail=F)

    int_power = stats::pchisq(q.thres,df=1,ncp=stats::qchisq(int_pval,df=1,lower.tail=F),lower.tail=F)
    a_power = stats::pchisq(q.thres,df=1,ncp=stats::qchisq(a_pval,df=1,lower.tail=F),lower.tail=F)
    b_power = stats::pchisq(q.thres,df=1,ncp=stats::qchisq(b_pval,df=1,lower.tail=F),lower.tail=F)
  } else
  {
    int_power=NA; a_power=NA; b_power=NA
  }
  return(c("power_int"=int_power,"power_snp1"=a_power,"power_snp2"=b_power))
}

