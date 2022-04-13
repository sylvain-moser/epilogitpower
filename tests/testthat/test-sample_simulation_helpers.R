context("Unit_test_for functions generate the simulation samples")
test_that("make_grr_correct", {
  actual=make_GRR_matrix(a1 =1,a2 =1,b1=1,b2=1,i=1.5)
  expected=matrix(c(1,1,1,1,1.5,2.25,1,2.25,5.0625),ncol=3) # from the epipower website
  expect_equal(actual,expected)
})


test_that("compute_gen_frequencies",{
  actual=compute_geno_freqs(mafa = 0.1,mafb = 0.1)
  expected=matrix(c(0.6561,0.1458,0.0081,0.1458,0.0324,0.0019,0.0081,0.0018,0.0001),ncol=3)
  expect_equal(actual,expected,tolerance=0.01)
})


test_that("compute_ref_geno_penetrance",{
  actual=compute_ref_geno_penetrance(prevalence = 0.1,grr = matrix(c(1,1,1,1,1.5,2.25,1,2.25,5.0625),ncol=3),geno_freqs = matrix(c(0.6561,0.1458,0.0081,0.1458,0.0324,0.0019,0.0081,0.0018,0.0001),ncol=3))
  expected=0.097933
  expect_equal(actual,expected,tolerance=0.01)
})


test_that("computed_ref_geno_penetrance_equald_inputed",{
  pD_reg_geno=0.05
  GRRs=ORs_to_GRRs(make_ORs_matrix(1.1,1.2,1.05,"dosage"),penetrance_ref_genotype = pD_reg_geno)
  pD=sum(GRRs*pD_reg_geno*matrix(c(0.6561,0.1458,0.0081,0.1458,0.0324,0.0019,0.0081,0.0018,0.0001),ncol=3))
  actual=compute_ref_geno_penetrance(prevalence = pD,grr = GRRs,geno_freqs = matrix(c(0.6561,0.1458,0.0081,0.1458,0.0324,0.0019,0.0081,0.0018,0.0001),ncol=3))
  expect_equal(actual,pD_reg_geno,tolerance=0.01)
})



test_that("compute_geno_penetrace",{
  actual=compute_geno_penetrances(grr = matrix(c(1,1,1,1,1.5,2.25,1,2.25,5.0625),ncol=3),ref_geno_penetrance = 0.097933)
  expected=matrix(c(0.0979,0.0979,0.0979,0.0979,0.1469,0.2203,0.0979,0.2203,0.4958),ncol=3)
  expect_equal(actual,expected,tolerance=0.01)
})

test_that("compute_geno_ORs",{
  actual=compute_geno_ORs(geno_penetrances = matrix(c(0.0979,0.0979,0.0979,0.0979,0.1469,0.2203,0.0979,0.2203,0.4958),ncol=3),ref_geno_penetrance = 0.097933)
  expected=matrix(c(1,1,1,1,1.5861,2.6033,1,2.6033,9.0571),ncol = 3)
  expect_equal(actual,expected,tolerance=0.01)
  })



test_that("to_geno_ORS_works",{
  actual=to_geno_OR(prevalence = 0.05,mafa = 0.2,mafb = 0.15,grr = make_grr_from_arr(a=0.7,b=1.3,i=1.2))
  expected=matrix(c(1,0.6890,0.4770,1.3211,1.0974,0.9132,1.7546,1.7700,1.7855),ncol = 3)
  expect_equal(actual,expected,tolerance=0.001)
})


test_that("geno_to_df_produce_right_dosages",{
  actual=geno_to_dosage(data.frame("pheno"=c(1,1,1),geno=c("AABb","AaBB","AaBb")))
  expected=data.frame("pheno"=c(1,1,1),"snp1"=c(0,1,1),"snp2"=c(1,0,1))
  expect_equal(actual,expected)
})
