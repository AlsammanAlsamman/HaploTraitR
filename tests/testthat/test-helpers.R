test_that("genotype calls are standardised", {
  expect_equal(HaploTraitR:::normalize_calls(c("A", "R", "AG", "N", "--", "A/G", NA, "AT")),
               c("AA", "AG", "AG", "NN", "NN", "AG", "NN", "AT"))
})

test_that("allele dosages count the minor allele and keep missing calls", {
  g <- matrix(c("AA", "AG", "GG", "AA", "NN"), nrow = 1)
  expect_equal(as.vector(convertGenoBi2Numeric(g)), c(0L, 1L, 2L, 0L, NA))
})

test_that("LD r2 equals the squared dosage correlation", {
  d <- rbind(c(0, 0, 2, 2, 1), c(0, 0, 2, 2, 1), c(2, 0, 0, 2, 0))
  r2 <- HaploTraitR:::ld_r2(d)
  expect_equal(r2[1, 2], 1)
  expect_equal(r2[1, 3], cor(d[1, ], d[3, ])^2)
})

test_that("nearby SNPs are found within the threshold", {
  res <- find_nearby_snps(c(1000, 3000), c(1500, 2500, 5000), 1000)
  expect_equal(res, list(1500, 2500))
})

test_that("LD components group linked SNPs", {
  adj <- matrix(FALSE, 4, 4)
  adj[1, 2] <- adj[2, 1] <- TRUE
  adj[3, 4] <- adj[4, 3] <- TRUE
  expect_equal(HaploTraitR:::ld_components(adj), c(1L, 1L, 2L, 2L))
})

test_that("compact letter display separates different groups", {
  letters <- HaploTraitR:::cld_letters(c("H3", "H2", "H1"), cbind(c("H3", "H3"), c("H2", "H1")))
  expect_equal(unname(letters), c("a", "b", "b"))
  none <- HaploTraitR:::cld_letters(c("H1", "H2"), NULL)
  expect_equal(unname(none), c("a", "a"))
})

test_that("phenotypes are cleaned and replicates averaged", {
  p <- suppressMessages(prepare_pheno(data.frame(id = c("a", "a", "b", "c"), y = c(1, 3, NA, 5))))
  expect_equal(p$Sample, c("a", "c"))
  expect_equal(p$Pheno, c(2, 5))
})

test_that("configuration round-trips through a text file", {
  on.exit(reset_config())
  suppressMessages({
    set_config(list(ld_threshold = 0.5, higher_is_better = FALSE))
    f <- tempfile(fileext = ".txt")
    save_config(f)
    reset_config()
    load_config(f)
  })
  expect_equal(get_config("ld_threshold"), 0.5)
  expect_false(get_config("higher_is_better"))
  expect_null(get_config("phenotypeunit"))
  expect_warning(set_config(list(not_a_parameter = 1)))
})
