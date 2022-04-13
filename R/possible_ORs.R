
#' Prediciton of possible interaction OR.
#' @description  Predicts the range of possible ORs for the interaction given data about the marginal effect of both SNPs
#' and epidemiological data.
#'
#' @param penetrance_ref_genotype The prevalence of the reference genotype
#' @param prevalence_up The upper bound of the range of possible prevalences. Default=1
#' @param prevalence_low The lower bound of the range of possible prevalences. Default=0
#' @param int_OR_upper The upper bound of the range of possible ORs for the interaction. Default=2
#' @param int_OR_lower The lower bound of the range of possible ORs for the interaction. Default=0.5
#' @param mafa maf of the snp1
#' @param mafb maf of the snp2
#' @param a_OR Odds Ratio for the additive effect of snp1
#' @param b_OR Odds Ratio for the additive effect of snp2
#' @param encoding Genetic encoding models for the two snps
#' @param plot Wheter to return a plot of the predicted prevalence in function of the interaction OR.
#'
#' @return a list of the maximum and minimum OR compatible with the upper and lower bound of the prevalence. If plot=T also returns a plot
#' of the predicted prevalence in function of the interaction OR.
#'
#' @export
#'
#' @examples
#' possible_int_OR(penetrance_ref_genotype = 0.0866,mafa = 0.22,mafb = 0.015,a_OR = 1.03,b_OR = 1.46,plot = FALSE,int_OR_upper = 2)
#' possible_int_OR(penetrance_ref_genotype = 0.0866,mafa = 0.22,mafb = 0.015,a_OR = 1.03,b_OR = 1.46,plot = TRUE,int_OR_upper = 2)

possible_int_OR=function(penetrance_ref_genotype=NULL,prevalence_up=1,prevalence_low=0,int_OR_upper=2,int_OR_lower=0.5,mafa,mafb,a_OR,b_OR,encoding="dosage",plot=FALSE){
prevalences=sapply(seq(int_OR_lower,int_OR_upper,0.01),prevalence_from_int_OR,penetrance_ref_genotype=penetrance_ref_genotype,mafa=mafa,mafb=mafb,a_OR=a_OR,b_OR=b_OR,encoding=encoding,simplify = T)
names(prevalences)=seq(int_OR_lower,int_OR_upper,0.01)
prevalences=prevalences[which(!is.na(prevalences))]
possible_ORs=prevalences[(prevalences <= prevalence_up) & (prevalences >=prevalence_low)]
if (length(prevalences[(prevalences <= prevalence_up)]) !=0){
  max_OR=as.numeric(names(which(possible_ORs==max(possible_ORs))))
} else {
  stop("The specified parameters yield a prevalence higher than its specified upper limit for any interaction OR")
}
if (length(prevalences[(prevalences >=prevalence_low)]) !=0){
  min_OR=as.numeric(names(which(possible_ORs==min(possible_ORs))))
} else{
  stop("The specified parameters yield a prevalence lower than its specified lower  limit for any interaction OR")
}
if (plot==TRUE){
  OR_plot=plot(x=names(possible_ORs),y=possible_ORs,main="Possible Interaction ORs",ylab ="Disease prevalence",xlab="Interaction OR")
  return(list("max_OR"=max_OR,"min_OR"=min_OR))
} else {
    return(list("max_OR"=max_OR,"min_OR"=min_OR))
  }
}

prevalence_from_int_OR=function(int_OR,penetrance_ref_genotype,mafa,mafb,a_OR,b_OR,encoding){
  ORs_matrix=make_ORs_matrix(a_OR=a_OR,b_OR=b_OR,int_OR = int_OR,encoding=encoding)
  grr_matrix=ORs_to_GRRs(ORs_matrix = ORs_matrix,penetrance_ref_genotype = penetrance_ref_genotype)
  tryCatch({
    prevalence=compute_prevalence(penetrance_ref_genotype,mafa,mafb,grr_matrix)
    return(prevalence)
  }, error=function(err){
    return(NA)
  })
}


compute_prevalence=function(penetrance_ref_genotype,mafa,mafb,grr){
  geno_freqs=compute_geno_freqs(mafa = mafa,mafb = mafb)
  geno_penetrances=compute_geno_penetrances(grr = grr,ref_geno_penetrance = penetrance_ref_genotype)
  prevalence=sum(geno_penetrances*geno_freqs)
  return(prevalence)
}
