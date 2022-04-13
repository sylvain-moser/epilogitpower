test_that("geno_probas_from_refGenoPenetrance match actual values",{
  actual=geno_probas_from_refGenoPenetrance(penetrance_ref_genotype=0.09496657,mafa = 0.1,mafb = 0.1,grr=matrix(c(1,1.15,1.09,1.1,1.4168,1.5040,1.05,1.5147,1.8009),ncol=3,byrow = T))
  expected=list("cases_geno_probs"=matrix(c(0.6231,0.1523 ,0.0081,0.1592 ,0.0436 ,0.0026,0.0084 ,0.0026 ,0.0002),ncol = 3),controls_geno_probs=matrix(c(0.6598,0.1451,0.0081,0.1443,0.0312,0.0017,0.0081,0.0017,0.001),ncol = 3))
  expect_equal(actual,expected,tolerance=0.0015)
})

test_that("geno_probas_from_OR_outputs_a_proba_matrix",{
  expect_type(geno_probas_from_OR(penetrance_ref_genotype=0.1,mafa = 0.1,mafb = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.04,encoding="dosage"),"list")
  cases_geno_probs=geno_probas_from_OR(penetrance_ref_genotype=0.7,mafa = 0.1,mafb = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.04,encoding="dosage")$cases_geno_probs
  controls_geno_probs=geno_probas_from_OR(penetrance_ref_genotype=0.7,mafa = 0.1,mafb = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.04,encoding="dosage")$controls_geno_probs
  expect_equal(class(cases_geno_probs),"matrix")
  expect_equal(class(controls_geno_probs),"matrix")
  expect_equal(sum(cases_geno_probs),1)
  expect_equal(sum(controls_geno_probs),1)
})

test_that("create_dataset_by_geno_sampling_from_OR_outputs_a_dosage_df",{
  dosage_df=create_dataset_by_geno_sampling_from_OR(n_cases = 1000,n_controls = 1000,penetrance_ref_genotype=0.09,mafa = 0.1,mafb = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,encoding = "dosage")
  expect_equal(dim(dosage_df),c(2000,3))
  expect_equal(length(dosage_df[dosage_df$pheno >1 | dosage_df$pheno< 0 ]),0)
  expect_equal(length(dosage_df[dosage_df$snp1 >2 | dosage_df$snp1 <0]),0)
  expect_equal(length(dosage_df[dosage_df$snp2 >2 | dosage_df$snp2 <0]),0)
})
