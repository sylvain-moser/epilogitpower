context("Unit_test_for functions performing Monte Carlo simualtions")
test_that("compute_pval_log_model_returns_a_pvalue",{
  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_snp1"],1)
  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR")[,"pvalue_snp1"],1)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_snp1"],0)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR")[,"pvalue_snp1"],0)

  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_snp2"],1)
  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR")[,"pvalue_snp2"],1)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_snp2"],0)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR")[,"pvalue_snp2"],0)

  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_int"],1)
  expect_lte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR")[,"pvalue_int"],1)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR")[,"pvalue_int"],0)
  expect_gte(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR")[,"pvalue_int"],0)

})

test_that("compute_pval_log_model_returns_error_when_arg_missing",{
  expect_error(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",input_type="GRR"))
  expect_error(compute_pval_log_model(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,prevalence = 0.1,encoding="dosage",input_type="OR"))
})

test_that("MC_power_analysis_return_a_power",{
  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_snp1"],1)
  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype=0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_snp1"],1)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_snp1"],0)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_snp1"],0)

  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_snp2"],1)
  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_snp2"],1)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_snp2"],0)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_snp2"],0)

  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_int"],1)
  expect_lte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_int"],1)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a1=1.2,a2=0.9,b1=1.14,b2=1.02,i=1.23,prevalence = 0.1,encoding="dosage",input_type="GRR",alpha = 0.05,n_simulations = 10)["power_int"],0)
  expect_gte(MC_power_analysis(n_controls=1000,n_cases=1000,maf_snp1 = 0.1,maf_snp2 = 0.1,a_OR=0.9,b_OR=1.1,int_OR = 1.5,prevalence = 0.1,encoding="dosage",penetrance_ref_genotype = 0.09,input_type="OR",alpha = 0.05,n_simulations = 10)["power_int"],0)
})


