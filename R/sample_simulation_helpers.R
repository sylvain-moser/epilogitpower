
### 2. Allele and genotype frequencies
### NOTE: the alleles that increases disease risk are always allele 'a' (locus A) and 'b' (locus B)
### Alleles 'a' and 'b' have frequencies pA and pB. Alleles 'A' and 'B' have frequencies qA and qB.
###
### 3. Two-locus genotype frequencies (NOTE: Locus A as rows and Locus B as columns)
### i.e. P(G)
###
###    | 	BB			Bb			bb		|
### ---+--------------------------------------------------+-----
### AA | 	f(AA).f(BB)		f(AA).f(Bb)		f(AA).f(bb)	|  f(AA)
### Aa | 	f(Aa).f(BB)		f(Aa).f(Bb)		f(Aa).f(bb)	|  f(Aa)
### aa | 	f(aa).f(BB)		f(aa).f(Bb)		f(aa).f(bb)	|  f(aa)
### ------------------------------------------------------------
###    |  f(BB)			f(Bb)			f(bb)		|  1

compute_geno_freqs=function(mafa,mafb){
  ### Locus A
  qA = 1 - mafa
  gen.freq_A = c()
  gen.freq_A[1] = qA ^ 2			## f(AA)
  gen.freq_A[2] = 2 * qA * mafa		## f(Aa)
  gen.freq_A[3] = mafa ^ 2			## f(aa)

  ### Locus B
  qB = 1 - mafb
  gen.freq_B = c()
  gen.freq_B[1] = qB ^ 2			## f(BB)
  gen.freq_B[2] = 2 * qB * mafb		## f(Bb)
  gen.freq_B[3] = mafb ^ 2			## f(bb)

  p.G = matrix(nrow=3,ncol=3)
  for (i in 1:3)
  {
    for (j in 1:3)
    {
      p.G[i,j] = gen.freq_A[i] * gen.freq_B[j]
    }
  }
  return (p.G)
}

###
### 4. Calculate penetrance of reference genotype
###    P(D|G=AABB) given P(D) and GRR table
###

###     P(D) = P(D|G=AABB).P(AABB) + ... + P(D|G=aabb).P(aabb)
### <=> P(D) / P(D|G=AABB) = P(AABB) + ... + GRR(aabb).P(aabb)
### <=> P(D|G=AABB) = P(D) / ( P(AABB) + ... + GRR(aabb).P(aabb) )

compute_ref_geno_penetrance=function(prevalence,grr,geno_freqs){
  return(prevalence/sum(grr*geno_freqs))
}


###
### 4a. Calculate genotype penetrances, P(D|G) from GRR table and P(D|G=AABB)
###    eg: GRR(AaBb) = P(D|G=AaBb)/P(D|G=AABB)
###    <=> P(D|G=AaBb) = GRR(AaBb) * P(D|G=AABB)
###

compute_geno_penetrances=function(grr,ref_geno_penetrance){
  p.d.G = grr * ref_geno_penetrance
  if (any(p.d.G > 1)) stop("Error: the GRR, allele frequencies and disease penetrance specified result in at least one penetrance > 1. Try to re-specify the model.")
  return(p.d.G)
}

###
### 5. Calculate genotype Odds ratios (OR)
###         OR(G) = [ P(D|G)/P(noD|G) ] / [ P(D|G=AABB)/P(noD|G=AABB) ]
###    <=>  OR(G) = [ P(D|G)/(1-P(D|G)) ] / [ P(D|G=AABB)/(1-P(D|G=AABB)) ]

compute_geno_ORs=function(geno_penetrances,ref_geno_penetrance){
  or.G = ( geno_penetrances / (1-geno_penetrances) ) / ( ref_geno_penetrance / (1-ref_geno_penetrance) )
  return(or.G)
}

###
### 6. Compute genotypes frequencies in the cases and controls:
### From P(D|G), P(G) and P(D) use Bayes' Theorem to get P(G|D) ( and similarly P(G|no D) )
###    P(G|D) will give me the expected frequencies of the 9 genotypes in cases, whereas P(G|no D)
###    are the expected frequencies in controls
###
###      P(G|D) = P(G).P(D|G) / P(D)
### and
###      P(G|no D) = P(G).P(no D|G) / P(no D)
### <=>  P(G|no D) = P(G).[1-P(D|G)] / [1-P(D)]

geno_probas=function(prevalence,geno_penetrances,geno_freqs){
  cases_geno_probs = (geno_freqs * geno_penetrances) / prevalence
  controls_geno_probs= (geno_freqs * (1-geno_penetrances)/ (1-prevalence))
  return (list("cases_geno_probs"=cases_geno_probs,"controls_geno_probs"=controls_geno_probs))
}




ORs_to_GRRs=function(ORs_matrix,penetrance_ref_genotype){
  apply(ORs_matrix,MARGIN = 1:2,function(OR) {OR/(1-penetrance_ref_genotype+(penetrance_ref_genotype*OR))})
}


make_ORs_matrix=function(a_OR,b_OR,int_OR,encoding){
  ors = matrix(ncol=3,nrow=3)

  ###
  ### 1. ORs (defined by a,b and i) for the dosage model
  ###
  ###    | 	BB			Bb			bb		|
  ### ---+--------------------------------------------------+-----
  ### AA | 	1			b			b^2		|  mean_AA
  ### Aa | 	a			a * b * i 		a*b^2*i^2	|  mean_Aa
  ### aa | 	a^2			a^2*b*i^2 		a^2*b^2*i^4	|  mean_aa
  ### ------------------------------------------------------------
  ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  if(encoding=="dosage"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = a_OR**2	# aaBB

    ors[1,2] = b_OR			# AABb
    ors[2,2] = a_OR * b_OR * int_OR		# AaBb
    ors[3,2] = a_OR**2 * b_OR *int_OR**2	# aaBb

    ors[1,3] = b_OR**2			# AAbb
    ors[2,3] = a_OR * b_OR**2 * int_OR**2	# Aabb
    ors[3,3] = a_OR**2 * b_OR**2 * int_OR**4	# aabb

    ###
    ### 2. ORs (defined by a,b and i) for the dominant model
    ###
    ###    | 	BB			Bb			bb		|
    ### ---+--------------------------------------------------+-----
    ### AA | 	1			b			    b		|  mean_AA
    ### Aa | 	a			a*b*i 		a*b*i	|  mean_Aa
    ### aa | 	a			a*b*i 		a*b*i	|  mean_aa
    ### ------------------------------------------------------------
    ###    |  mean_BB		mean_Bb		mean_bb	|  mean


  } else if (encoding=="dominant"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = a_OR	# aaBB

    ors[1,2] = b_OR			# AABb
    ors[2,2] = a_OR * b_OR * int_OR		# AaBb
    ors[3,2] = a_OR * b_OR *int_OR	# aaBb

    ors[1,3] = b_OR			# AAbb
    ors[2,3] = a_OR * b_OR * int_OR	# Aabb
    ors[3,3] = a_OR * b_OR * int_OR	# aabb

    ###
    ### 3. ORs (defined by a,b and i) for the recessive model
    ###
    ###    | 	BB			Bb			bb		|
    ### ---+--------------------------------------------------+-----
    ### AA | 	1			1			b		|  mean_AA
    ### Aa | 	1			1 		b	|  mean_Aa
    ### aa | 	a			a 		a*b*i	|  mean_aa ###
    ### ------------------------------------------------------------
    ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  } else if (encoding=="recessive"){
    ors[1,1] = 1			# AABB
    ors[2,1] = 1		# AaBB
    ors[3,1] = a_OR	# aaBB

    ors[1,2] = 1			# AABb
    ors[2,2] = 1		# AaBb
    ors[3,2] = a_OR	# aaBb

    ors[1,3] = b_OR			# AAbb
    ors[2,3] = b_OR	# Aabb
    ors[3,3] = a_OR * b_OR * int_OR	# aabb

    ###
    ### 4. ORs (defined by a,b and i) for the heterozygous model
    ###
    ###    | 	BB			Bb			bb		|
    ### ---+--------------------------------------------------+-----
    ### AA | 	1			b		  1		|  mean_AA
    ### Aa | 	a			a*b*i 	a	|  mean_Aa
    ### aa | 	1			b 		1	|  mean_aa ###
    ### ------------------------------------------------------------
    ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  } else if (encoding=="hetero"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = 1	# aaBB

    ors[1,2] = b_OR			# AABb
    ors[2,2] = a_OR * b_OR *int_OR		# AaBb
    ors[3,2] = b_OR	# aaBb

    ors[1,3] = 1			# AAbb
    ors[2,3] = a_OR	# Aabb
    ors[3,3] = 1	# aabb
  }

  ###
  ### 1. ORs (defined by a,b and i) for the dosage_recessive model
  ###
  ###    | 	BB			Bb			bb		|
  ### ---+--------------------------------------------------+-----
  ### AA | 	1			1			b		|  mean_AA
  ### Aa | 	a			a 		a*b*i	|  mean_Aa
  ### aa | 	a^2			a^2 		a^2*b*i^2 |  mean_aa   ## a=2, b=1 therefore the interaction term=2*1=2. There is two interaction terms and because of the multiplicative model (OR not logOR) OR_int needs to be squared
  ### ------------------------------------------------------------
  ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  else if (encoding=="dosage_recessive"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = a_OR**2	# aaBB

    ors[1,2] = 1			# AABb
    ors[2,2] = a_OR 		# AaBb
    ors[3,2] = a_OR**2	# aaBb

    ors[1,3] = b_OR			# AAbb
    ors[2,3] = a_OR*b_OR*int_OR	# Aabb
    ors[3,3] = a_OR**2 * b_OR * int_OR**2	# aabb
  }

  ###
  ### 6. ORs (defined by a,b and i) for the dominant (snp a)  recessive (snp b) model
  ###
  ###    | 	BB			Bb			bb		|
  ### ---+--------------------------------------------------+-----
  ### AA | 	1			1			b		|  mean_AA
  ### Aa | 	a			a 		a*b*i	|  mean_Aa
  ### aa | 	a		  a 		a*b*i |  mean_aa   ## not sure if need to square i
  ### ------------------------------------------------------------
  ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  else if (encoding=="dominant_recessive"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = a_OR	# aaBB

    ors[1,2] = 1			# AABb
    ors[2,2] = a_OR 		# AaBb
    ors[3,2] = a_OR	# aaBb

    ors[1,3] = b_OR			# AAbb
    ors[2,3] = a_OR*b_OR*int_OR	# Aabb
    ors[3,3] = a_OR* b_OR * int_OR	# aabb
  }
  ###
  ### 7. ORs (defined by a,b and i) for the het (snp a)  recessive (snp b) model
  ###
  ###    | 	BB			Bb			bb		|
  ### ---+--------------------------------------------------+-----
  ### AA | 	1			1			b		|  mean_AA
  ### Aa | 	a			a 		a*b*i	|  mean_Aa  ## square i ??
  ### aa | 	1		  1 		b |  mean_aa
  ### ------------------------------------------------------------
  ###    |  mean_BB		mean_Bb		mean_bb	|  mean

  else if (encoding=="hetero_recessive"){
    ors[1,1] = 1			# AABB
    ors[2,1] = a_OR		# AaBB
    ors[3,1] = 1	# aaBB

    ors[1,2] = 1			# AABb
    ors[2,2] = a_OR 		# AaBb
    ors[3,2] = 1	# aaBb

    ors[1,3] = b_OR			# AAbb
    ors[2,3] = a_OR*b_OR*int_OR	# Aabb
    ors[3,3] = b_OR	# aabb
  }

    return (ors)
}



#' @importFrom dplyr "%>%"
#' @noRD
geno_to_dosage=function(geno_df){
  df= geno_df %>%
      dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 2 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 2 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 2)) %>%
      dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 1 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 1 , geno=="AAbb" ~ 2, geno=="Aabb" ~ 2 ,geno=="aabb" ~ 2)) %>%
      dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}


#' @importFrom dplyr "%>%"
#' @noRD
geno_to_dominant=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 1 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 1 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 1 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 1 , geno=="AAbb" ~ 1, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#' @importFrom dplyr "%>%"
#' @noRD
geno_to_recessive=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 0, geno=="aaBB" ~ 1 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 0,geno=="aaBb" ~ 1 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 0 ,geno=="aabb" ~ 1)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 0,geno=="aaBb" ~ 0 , geno=="AAbb" ~ 1, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#' @importFrom dplyr "%>%"
#' @noRD
geno_to_hetero=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 0 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 0)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 1 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 1 , geno=="AAbb" ~ 0, geno=="Aabb" ~ 0 ,geno=="aabb" ~ 0)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#' @importFrom dplyr "%>%"
#' @noRD
geno_to_dosage_recessive=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 2 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 2 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 2)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 0,geno=="aaBb" ~ 0 , geno=="AAbb" ~ 1, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#' @importFrom dplyr "%>%"
#' @noRD
geno_to_dominant_recessive=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 1 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 1 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 0,geno=="aaBb" ~ 0 , geno=="AAbb" ~ 1, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#' @importFrom dplyr "%>%"
#' @noRD
geno_to_hetero_recessive=function(geno_df){
  df= geno_df %>%
    dplyr::mutate(snp1=dplyr::case_when(geno=="AABB" ~ 0, geno=="AaBB"~ 1, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 1,geno=="aaBb" ~ 0 ,geno== "AAbb" ~ 0, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 0)) %>%
    dplyr::mutate(snp2=dplyr::case_when(geno=="AABB" ~ 0,geno=="AaBB"~ 0, geno=="aaBB" ~ 0 , geno=="AABb" ~ 0 ,geno=="AaBb" ~ 0,geno=="aaBb" ~ 0 , geno=="AAbb" ~ 1, geno=="Aabb" ~ 1 ,geno=="aabb" ~ 1)) %>%
    dplyr::select(c("pheno","snp1","snp2"))
  return(df)
}

#### test helper function: these function allows to use epipower results to test my functions:
###
### 1. GRRs (defined by a,b and i)
###
###    | 	BB			Bb			bb		|
### ---+--------------------------------------------------+-----
### AA | 	1			b			b^2		|  mean_AA
### Aa | 	a			a * b * i 		a*b^2*i^2	|  mean_Aa
### aa | 	a^2			a^2*b*i^2 		a^2*b^2*i^4	|  mean_aa
### ------------------------------------------------------------
###    |  mean_BB		mean_Bb		mean_bb	|  mean




make_grr_from_arr=function(a,b,i){
  grr = matrix(ncol=3,nrow=3)

  grr[1,1] = 1			# AABB
  grr[2,1] = a		# AaBB
  grr[3,1] = a**2			# aaBB

  grr[1,2] = b		# AABb
  grr[2,2] = a * b * i		# AaBb
  grr[3,2] = a**2 * b * i**2	# aaBb

  grr[1,3] = b**2		# AAbb
  grr[2,3] = a * b**2 * i**2	# Aabb
  grr[3,3] = a**2 * b**2 * i**4	# aabb

  return (grr)
}

make_GRR_matrix=function(a1,a2,b1,b2,i){
  grr = matrix(ncol=3,nrow=3)

  grr[1,1] = 1			# AABB
  grr[2,1] = a1			# AaBB
  grr[3,1] = a2			# aaBB

  grr[1,2] = b1			# AABb
  grr[2,2] = a1 * b1 * i		# AaBb
  grr[3,2] = a2 * b1 * i^2	# aaBb

  grr[1,3] = b2			# AAbb
  grr[2,3] = a1 * b2 * i^2	# Aabb
  grr[3,3] = a2 * b2 * i^4	# aabb

  return (grr)
}




#This function is here only to check the transformation of OR to GRR
to_geno_OR=function(prevalence,mafa,mafb,grr){
  geno_freqs=compute_geno_freqs(mafa = mafa,mafb = mafb)
  ref_geno_penetrance=compute_ref_geno_penetrance(prevalence = prevalence,grr = grr,geno_freqs = geno_freqs)
  geno_penetrances=compute_geno_penetrances(grr = grr,ref_geno_penetrance = ref_geno_penetrance)
  geno_ORs=compute_geno_ORs(geno_penetrances = geno_penetrances,ref_geno_penetrance = ref_geno_penetrance)
  return(geno_ORs)
}
