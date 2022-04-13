test_that("geno probas match actual values",{
  actual=geno_probas(prevalence = 0.1,geno_penetrances = matrix(c(0.0979,0.0979,0.0979,0.0979,0.1469,0.2203,0.0979,0.2203,0.4958),ncol=3),geno_freqs = matrix(c(0.6561,0.1458,0.0081,0.1458,0.0324,0.0019,0.0081,0.0018,0.0001),ncol=3) )
  expected=list("cases_geno_probs"=matrix(c(0.6425,0.1428,0.0079,0.1428,0.0476,0.0040,0.0079,0.0040,0.0005),ncol = 3),controls_geno_probs=matrix(c(0.6576,0.1477,0.0083,0.1445,0.0308,0.0016,0.0079,0.0015,0.001),ncol = 3))
  expect_equal(actual,expected,tolerance=0.01)
})


test_that("geno_probas_from_prevalence match actual values",{
  actual=geno_probas_from_prevalence (prevalence = 0.1,mafa = 0.1,mafb = 0.1,grr=make_grr_from_arr(a =0.9,b=1.1,i=1.5))
  expected=list("cases_geno_probs"=matrix(c(0.6428,0.1286,0.0064,0.1571,0.0471,0.0035,0.0096,0.0043,0.0005),ncol = 3),controls_geno_probs=matrix(c(0.6576,0.1477,0.0083,0.1445,0.0308,0.0016,0.0079,0.0015,0.001),ncol = 3))
  expect_equal(actual,expected,tolerance=0.0015)
})

test_that("geno_probas outputs_a_proba_matrix",{
  expect_type(geno_probas_from_GRR(prevalence = 0.1,mafa = 0.1,mafb = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23),"list")
  cases_geno_probs=geno_probas_from_GRR(prevalence = 0.1,mafa = 0.1,mafb = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23)$cases_geno_probs
  controls_geno_probs=geno_probas_from_GRR(prevalence = 0.1,mafa = 0.1,mafb = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23)$controls_geno_probs
  expect_equal(class(cases_geno_probs),"matrix")
  expect_equal(class(controls_geno_probs),"matrix")
  expect_equal(sum(cases_geno_probs),1)
  expect_equal(sum(controls_geno_probs),1)
})


test_that("create_dataset_by_geno_sampling_from_GRR_outputs_a_dosage_df",{
  dosage_df=create_dataset_by_geno_sampling_from_GRR(n_cases = 1000,n_controls = 1000,prevalence = 0.1,mafa = 0.1,mafb = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,encoding = "dosage")
  expect_equal(dim(dosage_df),c(2000,3))
  expect_equal(length(dosage_df[dosage_df$pheno >1 | dosage_df$pheno< 0 ]),0)
  expect_equal(length(dosage_df[dosage_df$snp1 >2 | dosage_df$snp1 <0]),0)
  expect_equal(length(dosage_df[dosage_df$snp2 >2 | dosage_df$snp2 <0]),0)
})


