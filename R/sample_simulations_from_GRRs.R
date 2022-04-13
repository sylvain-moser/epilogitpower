geno_probas_from_GRR=function(prevalence,mafa,mafb,a1,a2,b1,b2,i){
  grr_matrix=make_GRR_matrix(a1,a2,b1,b2,i)
  geno_probas=geno_probas_from_prevalence(prevalence=prevalence,mafa = mafa,mafb = mafb,grr = grr_matrix)
  return(geno_probas)
}

geno_probas_from_prevalence=function(prevalence,mafa,mafb,grr){
  geno_freqs=compute_geno_freqs(mafa = mafa,mafb = mafb)
  ref_geno_penetrance=compute_ref_geno_penetrance(prevalence = prevalence,grr = grr,geno_freqs = geno_freqs)
  geno_penetrances=compute_geno_penetrances(grr = grr,ref_geno_penetrance = ref_geno_penetrance)
  assertthat::assert_that(all(all(geno_penetrances*geno_freqs >=0) & all(geno_penetrances*geno_freqs <=1) & (sum(geno_penetrances*geno_freqs) <=1) & (sum(geno_penetrances*geno_freqs) >=0)),msg="Error: the GRRs, allele frequencies and prevalence specified resulted in at least one  genotype penetrance > 1 or the prevalence >1. Try to re-specify the model.")
  stratified_geno_probas=geno_probas(prevalence = prevalence,geno_penetrances = geno_penetrances,geno_freqs = geno_freqs)
  return(stratified_geno_probas)
}


create_dataset_by_geno_sampling_from_GRR=function(n_cases,n_controls,prevalence,mafa,mafb,a1,a2,b1,b2,i,encoding){
  stratified_geno_proba=geno_probas_from_GRR(prevalence=prevalence,mafa=mafa,mafb=mafb,a1=a1,a2=a2,b1=b1,b2=b2,i=i)
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

