test_that("irregular CpGs map exactly and drive true-locus aggregation", {
  raw_reads <- data.table::data.table(
    chr = "chr1", read_start = 100, read_end = 220,
    original_hap_string = "010", count = 2, strand = "+",
    sample = "sample1", annotation_key = "annotation1"
  )
  cache <- list(annotation1 = list(chr1 = c(100, 135, 220, 500)))

  mapped <- MscoreDMR:::.map_reads_to_cpg(raw_reads, cache)
  loci <- MscoreDMR:::.build_master_loci(mapped$reads, cache)

  expect_true(mapped$reads$mapping_valid)
  expect_equal(mapped$reads$cpg_index_start, 1L)
  expect_equal(mapped$reads$cpg_index_end, 3L)
  expect_equal(loci$pos, c(100, 135, 220))
  expect_equal(loci$T, c(2, 2, 2))
  expect_equal(loci$S, c(2, 2, 2))
  expect_false(500 %in% loci$pos)
})

test_that("regional weights use local CpGs and the full-read global state", {
  raw_reads <- data.table::data.table(
    chr = "chr1", read_start = 100, read_end = 300,
    original_hap_string = "100", count = 5, strand = "+",
    sample = "sample1", annotation_key = "annotation1"
  )
  cache <- list(annotation1 = list(chr1 = c(100, 200, 300)))
  reads <- MscoreDMR:::.map_reads_to_cpg(raw_reads, cache)$reads
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(200, 300))

  matrices <- MscoreDMR:::.build_region_matrices(
    reads, cache, region, "sample1"
  )

  expect_true(reads$is_methylated_global)
  expect_equal(unname(matrices$T_matrix[1, 1]), 10)
  expect_equal(unname(matrices$S_matrix[1, 1]), 10)
  expect_equal(unname(matrices$N_sq_matrix[1, 1]), 20)
})

test_that("inclusive boundaries are exact and interval overlap is insufficient", {
  raw_reads <- data.table::data.table(
    chr = "chr1", read_start = 100, read_end = 300,
    original_hap_string = "11", count = 1, strand = "+",
    sample = "sample1", annotation_key = "annotation1"
  )
  cache <- list(annotation1 = list(chr1 = c(100, 300)))
  reads <- MscoreDMR:::.map_reads_to_cpg(raw_reads, cache)$reads

  boundary_region <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100, 300)
  )
  boundary <- MscoreDMR:::.build_region_matrices(
    reads, cache, boundary_region, "sample1"
  )
  expect_equal(unname(boundary$T_matrix[1, 1]), 2)
  expect_equal(unname(boundary$N_sq_matrix[1, 1]), 4)

  empty_region <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(150, 250)
  )
  empty <- MscoreDMR:::.build_region_matrices(
    reads, cache, empty_region, "sample1"
  )
  expect_length(empty$regions, 0L)
  expect_equal(empty$qc$reads_without_local_cpg, 1L)
})

test_that("invalid strings and CpG length mismatches are discarded with QC", {
  raw_reads <- data.table::data.table(
    chr = c("chr1", "chr1"),
    read_start = c(100, 100), read_end = c(300, 300),
    original_hap_string = c("01", "0X1"),
    count = c(1, 1), strand = c("+", "+"),
    sample = c("sample1", "sample1"),
    annotation_key = c("annotation1", "annotation1")
  )
  cache <- list(annotation1 = list(chr1 = c(100, 200, 300)))

  mapped <- NULL
  expect_message(
    mapped <- MscoreDMR:::.map_reads_to_cpg(raw_reads, cache),
    "Discarded 2"
  )

  expect_equal(nrow(mapped$reads), 0L)
  expect_equal(mapped$qc$cpg_length_mismatch_records, 1L)
  expect_equal(mapped$qc$invalid_haplotype_string_records, 1L)
  expect_equal(mapped$qc$invalid_mhap_records, 2L)
})

test_that("CpG files support shared paths and strict named sample matching", {
  cpg1 <- tempfile(fileext = ".tsv")
  cpg2 <- tempfile(fileext = ".tsv")
  writeLines("chr1\t100", cpg1)
  writeLines("chr1\t200", cpg2)
  on.exit(unlink(c(cpg1, cpg2)), add = TRUE)

  shared <- MscoreDMR:::.normalise_cpg_files(cpg1, c("s1", "s2"))
  expect_equal(names(shared), c("s1", "s2"))
  expect_equal(shared[1], shared[2], ignore_attr = TRUE)

  per_sample <- MscoreDMR:::.normalise_cpg_files(
    c(s2 = cpg2, s1 = cpg1), c("s1", "s2")
  )
  expect_match(unname(per_sample[1]), basename(cpg1), fixed = TRUE)
  expect_match(unname(per_sample[2]), basename(cpg2), fixed = TRUE)
  expect_error(
    MscoreDMR:::.normalise_cpg_files(c(cpg1, cpg2), c("s1", "s2")),
    "must have unique"
  )
})

test_that("multiple reads aggregate correctly for both single-locus and regional matrices", {
  raw_reads <- data.table::data.table(
    chr = rep("chr1", 2),
    read_start = c(100, 200), 
    read_end = c(300, 400),
    original_hap_string = c("100", "000"), 
    count = c(5, 10),                      
    strand = rep("+", 2),
    sample = rep("sample1", 2),
    annotation_key = rep("annotation1", 2)
  )
  cache <- list(annotation1 = list(chr1 = c(100, 200, 300, 400)))

  mapped <- MscoreDMR:::.map_reads_to_cpg(raw_reads, cache)
  reads <- mapped$reads

  loci <- MscoreDMR:::.build_master_loci(reads, cache)
  expect_equal(loci$pos, c(100, 200, 300, 400))

  expect_equal(loci[pos == 100]$T, 5)
  expect_equal(loci[pos == 100]$S, 5)
  expect_equal(loci[pos == 200]$T, 15)
  expect_equal(loci[pos == 200]$S, 5)
  expect_equal(loci[pos == 300]$T, 15)
  expect_equal(loci[pos == 300]$S, 5)
  expect_equal(loci[pos == 400]$T, 10)
  expect_equal(loci[pos == 400]$S, 0)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(200, 300))
  
  matrices <- MscoreDMR:::.build_region_matrices(
    reads, cache, region, "sample1"
  )

  expect_equal(unname(matrices$T_matrix[1, 1]), 10 + 20)      # 30
  expect_equal(unname(matrices$S_matrix[1, 1]), 10 + 0)       # 10
  expect_equal(unname(matrices$N_sq_matrix[1, 1]), 20 + 40)   # 60
})
