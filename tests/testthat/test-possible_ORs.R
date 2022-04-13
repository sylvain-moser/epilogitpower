test_that("possible_OR_outputs_two_floats",{

  expect_named(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5),c("max_OR","min_OR"))

  expect_type(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$max_OR,"double")
  expect_type(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$min_OR,"double")

  expect_lte(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$max_OR,10)
  expect_gte(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$max_OR,-10)
  expect_lte(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$min_OR,10)
  expect_gte(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=1.2,b_OR=1.4,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5)$min_OR,-10)

})

test_that("possible_OR_outputs_error_when_model_impossible",{
  expect_error(possible_int_OR(penetrance_ref_genotype=0.1,mafa = 0.8,mafb = 0.7,a_OR=4,b_OR=5,encoding="dosage",prevalence_low = 0,plot = F,prevalence_up =0.5))
})
