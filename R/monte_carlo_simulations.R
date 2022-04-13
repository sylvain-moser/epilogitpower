#compute_pval_log_model based on Wald-test of the coefficients
compute_pval_log_model=function(dummy_iterator,n_cases,n_controls,maf_snp1,maf_snp2,a_OR=NULL,b_OR=NULL,int_OR=NULL,a1=NULL,a2=NULL,b1=NULL,b2=NULL,i=NULL,prevalence=NULL,penetrance_ref_genotype=NULL,encoding="dosage",input_type){
  if (input_type=="OR"){
    #user input check:
    assertthat::assert_that(!is.null(n_cases) && !is.null(n_controls) && !is.null(maf_snp1) && !is.null(maf_snp2) && !is.null(a_OR) && !is.null(b_OR) && !is.null(int_OR) && !is.null(encoding) && !is.null(penetrance_ref_genotype),msg = "If input_type == 'OR' the following arguments have to be provided : n_cases,n_controls,maf_snp1,maf_snp2,a_OR,b_OR,int_OR,penetrance_ref_genotype,encoding")
    #
    df=create_dataset_by_geno_sampling_from_OR(n_cases=n_cases,n_controls=n_controls,mafa=maf_snp1,mafb=maf_snp2,a_OR=a_OR,b_OR=b_OR,int_OR=int_OR,penetrance_ref_genotype,encoding)
  } else if (input_type=="GRR"){
    #user input check:
    assertthat::assert_that(! is.null(n_cases) && !is.null(n_controls) && !is.null(maf_snp1) && !is.null(maf_snp2) && !is.null(a1) && !is.null(a2) && !is.null(b1) && !is.null(b2) && !is.null(i) && !is.null(prevalence) && !is.null(encoding),msg = "If input_type == 'GRR' the following arguments have to be provided : n_cases,n_controls,prevalence,maf_snp1,maf_snp2,a1,a2,b1,b2,i,encoding")
    #
    df=create_dataset_by_geno_sampling_from_GRR(n_cases=n_cases,n_controls=n_controls,prevalence=prevalence,mafa=maf_snp1,mafb=maf_snp2,a1=a1,a2=a2,b1=b1,b2=b2,i=i,encoding=encoding)
  }
  pval_df=tryCatch({
  mod=stats::glm(pheno ~ snp1*snp2,data=df,family = "binomial")
    pvalue_snp1=summary(mod)$coef["snp1","Pr(>|z|)"]
    pvalue_snp2=summary(mod)$coef["snp2","Pr(>|z|)"]
    pvalue_int=summary(mod)$coef["snp1:snp2","Pr(>|z|)"]
    return(data.frame("pvalue_snp1"=pvalue_snp1,"pvalue_snp2"=pvalue_snp2,"pvalue_int"=pvalue_int))
  },error=function(err){
    if (input_type=="OR"){ ### the trycatch system here aims to re-run the sampling (random) if the model fail to get other values. it does it twice now. Could implement to do this until succes?
      df=create_dataset_by_geno_sampling_from_OR(n_cases=n_cases,n_controls=n_controls,mafa=maf_snp1,mafb=maf_snp2,a_OR=a_OR,b_OR=b_OR,int_OR=int_OR,penetrance_ref_genotype,encoding)
    } else if (input_type=="GRR"){
      df=create_dataset_by_geno_sampling_from_GRR(n_cases=n_cases,n_controls=n_controls,prevalence=prevalence,mafa=maf_snp1,mafb=maf_snp2,a1=a1,a2=a2,b1=b1,b2=b2,i=i,encoding=encoding)
    }
    pval_df2=tryCatch({
    mod=stats::glm(pheno ~ snp1*snp2,data=df,family = "binomial")
    pvalue_snp1=summary(mod)$coef["snp1","Pr(>|z|)"]
    pvalue_snp2=summary(mod)$coef["snp2","Pr(>|z|)"]
    pvalue_int=summary(mod)$coef["snp1:snp2","Pr(>|z|)"]
    return(data.frame("pvalue_snp1"=pvalue_snp1,"pvalue_snp2"=pvalue_snp2,"pvalue_int"=pvalue_int))
    }, error=function(err){
      pvalue_snp1=NA
      pvalue_snp2=NA
      pvalue_int=NA
      return(data.frame("pvalue_snp1"=pvalue_snp1,"pvalue_snp2"=pvalue_snp2,"pvalue_int"=pvalue_int))
    }
    )
    return(pval_df2)
  }
  )
  return(pval_df)
}


#' Monte Carlo based power analysis
#'
#' Genotypes frequencies are computed separately for cases and
#' controls based on the prevalence of the diseases, the minor allele frequency of the 2 snps and their ORs.
#' n_cases cases and n_controls are then simulated using the genotype frequencies and used in a logistic regression
#' model to predict the diseases phenotype. The power is then computed by estimating the distribution of the Wald-test statistics under the alternative
#' by Monte-Carlo simulations.
#'
#'
#'
#' @param sample_size the total sample size if the sampling scheme is "naive"
#' @param n_cases the number of cases to sample if sampling scheme is "stratified"
#' @param n_controls the number of controls to sample if sampling scheme is "stratified"
#' @param maf_snp1 maf of the snp1
#' @param maf_snp2 maf of the snp2
#' @param a_OR Odds ratio for the additive effect of snp1
#' @param b_OR Odds ratio for the additive effect of snp2
#' @param int_OR Odds ratio for the interaction effect of snp1 and snp2
#' @param prevalence prevalence of the disease
#' @param sampling_method method used to generate the dataset. Either "naive" randomly drawing from a common population determined by the mafs or "Stratified" drawing cases and controls separately
#' @param alpha the alpha threshold for overall significance
#' @param n_simulations number of Monte Carlo simulations to run
#' @param threads number of threads to use for the MC simulations
#' @param encoding the genetic model encoding to choose from "dosage","dominant","recessive","hetero" with respect to minor allele. Only
#'
#' @return a sinlge float. The power
#' @export
#' @examples
#'MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)
#'MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR",alpha = 0.05,n_simulations = 10)

MC_power_analysis=function(n_cases,n_controls,maf_snp1,maf_snp2,a_OR=NULL,b_OR=NULL,int_OR=NULL,a1=NULL,a2=NULL,b1=NULL,b2=NULL,i=NULL,prevalence=NULL,penetrance_ref_genotype=NULL,encoding="dosage",input_type,alpha,n_simulations=100,threads=1){
  if (input_type=="OR"){
    #user input check:
    assertthat::assert_that(! is.null(n_cases) && !is.null(n_controls) && !is.null(maf_snp1) && !is.null(maf_snp2) && !is.null(a_OR) && !is.null(b_OR) && !is.null(int_OR) && !is.null(encoding) && !is.null(penetrance_ref_genotype),msg = "If input_type == 'OR' the following arguments have to be provided : n_cases,n_controls,maf_snp1,maf_snp2,a_OR,b_OR,int_OR,penetrance_ref_genotype,encoding")
  } else if (input_type=="GRR"){
    #user input check:
    assertthat::assert_that(! is.null(n_cases) && !is.null(n_controls) && !is.null(maf_snp1) && !is.null(maf_snp2) && !is.null(a1) && !is.null(a2) && !is.null(b1) && !is.null(b2) && !is.null(i) && !is.null(prevalence) && !is.null(encoding),msg = "If input_type == 'GRR' the following arguments have to be provided : n_cases,n_controls,prevalence,maf_snp1,maf_snp2,a1,a2,b1,b2,i,encoding")
  }
  pval_list=parallel::mclapply(seq(1:n_simulations),FUN = compute_pval_log_model,n_cases=n_cases,n_controls=n_controls,maf_snp1=maf_snp1,maf_snp2=maf_snp2,a_OR=a_OR,b_OR=b_OR,int_OR=int_OR,a1=a1,a2=a2,b1=b1,b2=b2,i=i,prevalence=prevalence,penetrance_ref_genotype=penetrance_ref_genotype,input_type=input_type,encoding=encoding,mc.cores = threads)
  pval_df=data.table::rbindlist(pval_list)
  pval_df=as.data.frame(pval_df)
  pval_df=pval_df[which(!is.na(pval_df$pvalue_snp1)),]
  power_snp1=length(pval_df$pvalue_snp1[pval_df$pvalue_snp1 <=alpha])/length(pval_df$pvalue_snp1)
  power_snp2=length(pval_df$pvalue_snp2[pval_df$pvalue_snp2 <=alpha])/length(pval_df$pvalue_snp2)
  power_int=length(pval_df$pvalue_int[pval_df$pvalue_int <=alpha])/length(pval_df$pvalue_int)
  return(c("power_snp1"=power_snp1,"power_snp2"=power_snp2,"power_int"=power_int))
}

