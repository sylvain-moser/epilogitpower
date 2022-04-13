test_that("compute_chi2_power_outputs_power",{
  expect_lte(compute_chi2_power(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05)["power_int"],1)
  expect_gte(compute_chi2_power(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05)["power_int"],0)

  expect_lte(compute_chi2_power(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.2,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR",alpha = 0.05)["power_snp1"],1)
  expect_gte(compute_chi2_power(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.2,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR",alpha = 0.05)["power_snp1"],0)

})

test_that("compute_pval_log_model_returns_error_when_arg_missing",{
  expect_error(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,a2,int_OR = 1.2,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR"))
  expect_error(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,encoding="dosage",input_type="GRR"))
})
