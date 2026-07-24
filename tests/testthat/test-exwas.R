test_that("exwas returns expected structure", {
  data(exposome_example)
  data(pheno_example)
  pheno_test <- pheno_example
  pheno_test$Disease_binary <- ifelse(pheno_test$Disease == "Yes", 1, 0)

  fit <- exwas(
    Disease_binary ~ omic_features + age + sex,
    pheno = pheno_test,
    omic_features = exposome_example,
    id_col = "exp_id",
    family = "binomial"
  )

  expect_type(fit, "list")
  expect_true("results" %in% names(fit))
  expect_true(all(c("variable","term","estimate","std.error","lcl","ucl","p.value") %in% names(fit$results)))
})
