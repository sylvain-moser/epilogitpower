### to estimate the genotype frequencies, starting from a pre-defined interaction OR:
geno_probas_from_OR=function(penetrance_ref_genotype,mafa,mafb,a_OR,b_OR,int_OR,encoding){
  ORs_matrix=make_ORs_matrix(a_OR=a_OR,b_OR=b_OR,int_OR = int_OR,encoding=encoding)
  grr_matrix=ORs_to_GRRs(ORs_matrix = ORs_matrix,penetrance_ref_genotype = penetrance_ref_genotype)
  geno_freqs=compute_geno_freqs(mafa = mafa,mafb = mafb)
  geno_probas=geno_probas_from_refGenoPenetrance(penetrance_ref_genotype = penetrance_ref_genotype,mafa=mafa,mafb = mafb,grr = grr_matrix)
  return(geno_probas)
}


geno_probas_from_refGenoPenetrance=function(penetrance_ref_genotype,mafa,mafb,grr){
  geno_freqs=compute_geno_freqs(mafa = mafa,mafb = mafb)
  geno_penetrances=compute_geno_penetrances(grr = grr,ref_geno_penetrance = penetrance_ref_genotype)
  prevalence=sum(geno_penetrances*geno_freqs)
  assertthat::assert_that(prevalence <=1 && prevalence >=0,msg="Error:the ORs, allele frequencies and reference allele penetrance specified resulted in the prevalence being > 1 or <0. Try to re-specify the model")
  assertthat::assert_that(all(all(geno_penetrances*geno_freqs >=0) & all(geno_penetrances*geno_freqs <=1) & (sum(geno_penetrances*geno_freqs) <=1) & (sum(geno_penetrances*geno_freqs) >=0)),msg="Error: the ORs, allele frequencies and reference allele penetrance specified resulted in at least one  genotype penetrance > 1 or the prevalence >1. Try to re-specify the model.")
  stratified_geno_probas=geno_probas(prevalence = prevalence,geno_penetrances = geno_penetrances,geno_freqs = geno_freqs)
  return(stratified_geno_probas)
}


create_dataset_by_geno_sampling_from_OR=function(n_cases,n_controls,penetrance_ref_genotype,mafa,mafb,a_OR,b_OR,int_OR,encoding){
  stratified_geno_proba=geno_probas_from_OR(penetrance_ref_genotype=penetrance_ref_genotype,mafa=mafa,mafb=mafb,a_OR=a_OR,b_OR=b_OR,int_OR=int_OR,encoding=encoding)
  cases_geno=data.frame("pheno"=rep(1,n_cases),geno=sample(x = c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb") ,size = n_cases,prob = as.vector(t(stratified_geno_proba$cases_geno_probs)),replace = T)) # maybe here should turn around!
  controls_geno=data.frame("pheno"=rep(0,n_controls),geno=sample(x = c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb") ,size = n_controls,prob = as.vector(t(stratified_geno_proba$controls_geno_probs)),replace = T))
  df=rbind(cases_geno,controls_geno)
  if (encoding=="geno"){
    return(df)
  } else if (encoding=="dosage"){
    df=geno_to_dosage(df)
    return(df)
  }else if (encoding=="dominant"){
    df=geno_to_dominant(df)
    return(df)
  }else if (encoding=="recessive"){
    df=geno_to_recessive(df)
    return(df)
  }else if (encoding=="hetero"){
    df=geno_to_hetero(df)
    return(df)
  }else if (encoding=="dosage_recessive"){
    df=geno_to_dosage_recessive(df)
    return(df)
  }else if (encoding=="dominant_recessive"){
    df=geno_to_dominant_recessive(df)
    return(df)
  }else if (encoding=="hetero_recessive"){
    df=geno_to_hetero_recessive(df)
    return(df)
  }
}



