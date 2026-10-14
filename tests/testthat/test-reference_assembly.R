test_that("reference_assembly curates lab reference strains and falls back to the organism", {
  expect_identical(reference_assembly(c("Bacillus subtilis", "Escherichia coli", "Homo sapiens")),
                   c("GCA_000009045", "GCA_000005845", "Homo sapiens"))
  expect_identical(formals(get_annotation)$assembly, quote(reference_assembly(organism)))
})
