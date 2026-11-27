test_that("plain Zymomonas mobilis resolves to the ZM4 reference, other organisms unchanged", {
  expect_identical(organism_alias(c("Zymomonas mobilis", " zymomonas  mobilis ", "Zymomonas mobilis ZM4", "Homo sapiens")),
                   c("Zymomonas mobilis ZM4", "Zymomonas mobilis ZM4", "Zymomonas mobilis ZM4", "Homo sapiens"))
  expect_identical(reference_folder_name("Zymomonas mobilis"), "zymomonas_mobilis_zm4")
  expect_identical(reference_folder_name("Zymomonas mobilis ZM4"), "zymomonas_mobilis_zm4")
  expect_identical(reference_folder_name("Homo sapiens"), "homo_sapiens")
  expect_identical(organism_alias(character()), character())
})
