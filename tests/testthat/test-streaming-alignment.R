test_that("difference events add repeated endpoints and cancellations exactly", {
  reads <- data.table::data.table(cpg_index_start = c(3L, 3L, 6L),
    cpg_index_end = c(6L, 5L, 6L), count = c(10, 5, 2^32),
    is_methylated_global = c(TRUE, FALSE, FALSE))
  delta <- list(S = numeric(7), T = numeric(7))
  both <- MscoreDMR:::.add_cpg_events(delta, reads)
  separate <- delta
  for (i in seq_len(nrow(reads))) separate <- MscoreDMR:::.add_cpg_events(separate, reads[i])
  expect_identical(both, separate)
  expect_equal(cumsum(both$S), c(0, 0, 10, 10, 10, 10, 0))
  expect_equal(cumsum(both$T), c(0, 0, 15, 15, 15, 10 + 2^32, 0))
})

test_that("streaming counts and indexed region matrices equal the old path", {
  root <- tempfile("stream-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  data.table::fwrite(data.frame(chr = rep(c("chr1", "chr2"), each = 4),
    pos = rep(c(100, 135, 220, 500), 2), end = rep(c(101, 136, 221, 501), 2)),
    cpg, sep = "\t", col.names = FALSE)
  # Unsorted records, identical records crossing batch boundaries, both strands,
  # last CpG, globally methylated but locally unmethylated, invalid data.
  a <- data.table::data.table(chr = c("chr2", rep("chr1", 6), "chrMissing"),
    start = c(100, 100, 500, 100, 135, 100, 100, 100),
    end = c(220, 220, 500, 220, 220, 220, 220, 100),
    hap = c("100", "100", "0", "100", "00", "01", "0X1", "1"),
    count = c(2, 5, 2^32, 5, 3, 1, 1, 1), strand = c("-", rep("+", 7)))
  b <- a[1:5]
  b <- b[!(chr == "chr1" & start == 500)]
  files <- stats::setNames(file.path(root, c("a.tsv", "b.tsv")), c("s1", "s2"))
  data.table::fwrite(a, files[1], sep = "\t", col.names = FALSE)
  data.table::fwrite(b, files[2], sep = "\t", col.names = FALSE)
  old <- MscoreDMR:::.prepare_mscore_data(files, cpg)
  regions <- GenomicRanges::GRanges(c("chr1", "chr1", "chr1", "chr2"),
    IRanges::IRanges(c(100, 136, 500, 135), c(500, 219, 500, 220)))
  expected <- MscoreDMR:::.build_region_matrices(old$master_read_dt, old$cpg_cache,
                                               regions, names(files))
  for (batch in c(1L, 2L, 100000L)) {
    store <- MscoreDMR:::.prepare_mscore_store(files, cpg, batch, root)
    actual <- MscoreDMR:::.store_locus_table(store)
    expect_equal(actual, old$master_locus_dt, tolerance = 0)
    expect_equal(store$qc, old$qc)
    got <- MscoreDMR:::.region_matrices_store(store, regions, names(files))
    expect_equal(got, expected, tolerance = 0)
    groups <- c(s1 = "A", s2 = "B")
    screened <- MscoreDMR:::.screen_mscore_store(store, groups, c("A", "B"), 30)
    expect_false(500 %in% screened$candidates[["chr1"]]$pos)
    expect_true(500 %in% MscoreDMR:::.store_locus_table(screened)$pos)
    # The excluded coordinate must STILL contribute to regional S/T/N_sq.
    expect_equal(MscoreDMR:::.region_matrices_store(screened, regions, names(files)), expected)
    unlink(store$root, recursive = TRUE)
  }
})

test_that("two-chromosome store smoothing shares exact global parameters", {
  root <- tempfile("smooth-store-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  loci <- data.table::CJ(chr = c("chr1", "chr2"), pos = seq(100, 2000, 100),
                         sample = paste0("s", 1:4))
  loci[chr == "chr2", pos := pos * 3]
  loci[, T := ifelse(chr == "chr1", 10, 200)]
  loci[, S := T * ifelse(sample %in% c("s1", "s2"), .8, .1)]
  loci[, M := S / T]
  groups <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")
  dense <- MscoreDMR:::.dense_loci_from_long(loci, names(groups))
  store <- MscoreDMR:::.screen_mscore_store(list(root = root, loci = dense, qc = list()),
                                          groups, c("A", "B"), 30)
  all_stats <- MscoreDMR:::.compute_candidate_smoothing_stats(loci, groups, c("A", "B"))
  expect_equal(store$cov_q75, attr(all_stats, "cov_q75"), tolerance = 0)
  expect_equal(store$density_span, attr(all_stats, "density_span"), tolerance = 0)
  old <- MscoreDMR:::.run_smoothing_pipeline(loci, groups, c("A", "B"), .1)
  new <- MscoreDMR:::.smooth_mscore_store(store, groups, c("A", "B"), .1,
                                        2500, 1000, 5, 1000, 30)
  expect_equal(as.data.frame(new), as.data.frame(old), tolerance = 0)
})

test_that("batch size is validated and private caches are cleaned on failure", {
  for (x in list(0, -1, NA_real_, Inf, 1.5, c(1, 2))) {
    expect_error(MscoreDMR:::.validate_chunk_nrows(x), "positive integer")
  }
  expect_identical(MscoreDMR:::.validate_chunk_nrows(100000), 100000L)
  root <- tempfile("bad-stream-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  bad <- file.path(root, "bad.tsv")
  writeLines("chr1\t100", cpg)
  writeLines("chr1\t100", bad)
  before <- list.files(root)
  expect_error(MscoreDMR:::.prepare_mscore_store(c(s1 = bad, s2 = bad), cpg,
                                               cache_dir = root), "six")
  expect_identical(list.files(root), before)
})

test_that("compressed inputs and per-sample annotations preserve identities", {
  root <- tempfile("per-annotation-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  annotations <- stats::setNames(file.path(root, c("a.cpg", "b.cpg")), c("a", "b"))
  writeLines(c("chr1\t100", "chr1\t150", "chr1\t200"), annotations[1])
  writeLines(c("chr1\t100", "chr1\t170", "chr1\t200"), annotations[2])
  files <- stats::setNames(file.path(root, c("a.mhap", "b.mhap")), c("a", "b"))
  for (path in files) writeLines(c("chr1\t100\t200\t100\t2\t-", "chr1\t200\t200\t0\t3\t+"), path)
  old <- MscoreDMR:::.prepare_mscore_data(files, annotations)
  packed <- stats::setNames(vapply(files, Rsamtools::bgzip, character(1)), names(files))
  store <- MscoreDMR:::.prepare_mscore_store(packed, annotations[c(2, 1)], 1L, root)
  empty <- MscoreDMR:::.region_matrices_store(store, GenomicRanges::GRanges(), names(files))
  expect_length(empty$regions, 0L)
  expect_equal(dim(empty$S_matrix), c(0L, 2L))
  actual <- MscoreDMR:::.store_locus_table(store)
  expect_equal(actual, old$master_locus_dt, tolerance = 0)
  regions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), c(200, 170)))
  expect_equal(MscoreDMR:::.region_matrices_store(store, regions, names(files)),
    MscoreDMR:::.build_region_matrices(old$master_read_dt, old$cpg_cache, regions, names(files)),
    tolerance = 0)
})

test_that("tabix shards keep large genomic coordinates in decimal notation", {
  root <- tempfile("large-position-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  mhap <- file.path(root, "reads.tsv")
  writeLines("chr1\t100000000", cpg)
  writeLines("chr1\t100000000\t100000000\t1\t10\t+", mhap)
  files <- c(s1 = mhap, s2 = mhap)
  store <- MscoreDMR:::.prepare_mscore_store(files, cpg, 1L, root)
  regions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100000000L, 100000000L))
  matrices <- MscoreDMR:::.region_matrices_store(store, regions, names(files))
  expect_equal(unname(matrices$S_matrix), matrix(10, 1, 2))
  expect_equal(unname(matrices$T_matrix), matrix(10, 1, 2))
  expect_equal(unname(matrices$N_sq_matrix), matrix(10, 1, 2))
})

test_that("the C reader reproduces the R reader exactly on sorted input", {
  root <- tempfile("c-reader-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  writeLines(paste0(rep(c("chr1", "chr2"), each = 6), "\t",
                    rep(c(100, 135, 220, 300, 301, 500), 2)), cpg)
  # Sorted by chromosome block and start, with ends out of order inside a
  # start, identical records split apart by another span, duplicate spans
  # with different strings or strands, every invalid-record category,
  # padded fields, a blank line and a CRLF line ending.
  lines <- c(
    "chr2\t100\t220\t101\t3\t+",
    "chr2\t100\t135\t10\t4\t-",
    "chr2\t100\t220\t101\t2\t+",
    "chr2\t100\t220\t001\t5\t+",
    "chr2\t100\t220\t101\t1\t-",
    "chr2\t135\t300\t0X0\t1\t+",
    "",
    "chr2\t135\t301\t0100\t2\t+",
    " chr2 \t220\t500\t1111\t6\t+\r",
    "chr1\t100\t100\t1\t2\t+",
    "chr1\t100\t300\t1\t2\t+",
    "chr1\t136\t219\t1\t2\t+",
    "chr1\t220\t500\t0000\t0\t+",
    "chr1\t220\t500\t0000\t7\t*",
    "chr1\t300\t220\t00\t7\t+",
    "chr1\t300\t301\tNA\t7\t+",
    "chr1\t301\t500\t01\t4294967296\t-",
    "chrMissing\t1\t2\t1\t1\t+")
  files <- stats::setNames(file.path(root, c("a.mhap", "b.mhap.gz")), c("a", "b"))
  writeLines(lines, files[["a"]])
  connection <- gzfile(files[["b"]], "w")
  writeLines(lines[-(1:5)], connection)
  close(connection)
  regions <- GenomicRanges::GRanges(c("chr1", "chr1", "chr2", "chr2", "chr2"),
    IRanges::IRanges(c(100, 136, 100, 136, 221), c(500, 300, 135, 400, 500)))
  r_store <- MscoreDMR:::.prepare_mscore_store(files, cpg, 2L, root,
                                               read_backend = "tabix")
  scan_calls <- 0L
  scanner <- getFromNamespace(".scan_mhap_file", "MscoreDMR")
  testthat::local_mocked_bindings(.scan_mhap_file = function(...) {
    scan_calls <<- scan_calls + 1L
    result <- scanner(...)
    expect_false(is.null(result))
    result
  }, .package = "MscoreDMR")
  c_store <- MscoreDMR:::.prepare_mscore_store(files, cpg, 2L, root)
  expect_equal(scan_calls, 2L)
  expect_equal(c_store$qc, r_store$qc, tolerance = 0)
  expect_equal(c_store$qc$total_mhap_records, 17 + 12)
  expect_equal(c_store$qc$invalid_haplotype_string_records, 2 * 2)
  expect_identical(c_store$loci, r_store$loci)
  expect_equal(MscoreDMR:::.region_matrices_store(c_store, regions, names(files)),
               MscoreDMR:::.region_matrices_store(r_store, regions, names(files)),
               tolerance = 0)
  # Without the read cache (candidate discovery only) the same loci result.
  counts_only <- MscoreDMR:::.prepare_mscore_store(files, cpg, 2L, root,
                                                   keep_reads = FALSE)
  expect_equal(counts_only$qc, r_store$qc, tolerance = 0)
  expect_length(counts_only$read_paths, 0L)
  for (store in list(r_store, c_store, counts_only)) unlink(store$root, recursive = TRUE)
})

test_that("unsorted input falls back to the R reader with identical results", {
  root <- tempfile("c-reader-unsorted-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  writeLines(paste0("chr1\t", c(100, 200, 300)), cpg)
  files <- stats::setNames(file.path(root, c("a.tsv", "b.tsv")), c("a", "b"))
  writeLines(c("chr1\t200\t300\t11\t2\t+", "chr1\t100\t200\t01\t3\t+"), files[1])
  writeLines(c("chr1\t100\t300\t011\t2\t+"), files[2])
  expect_message(store <- MscoreDMR:::.prepare_mscore_store(files, cpg, 1L, root),
                 "not sorted")
  expect_length(list.files(store$root, pattern = "^reads-1-"), 0L)
  reference <- MscoreDMR:::.prepare_mscore_store(files, cpg, 1L, root,
                                                 read_backend = "tabix")
  regions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 200), c(300, 300)))
  expect_equal(MscoreDMR:::.region_matrices_store(store, regions, names(files)),
               MscoreDMR:::.region_matrices_store(reference, regions, names(files)),
               tolerance = 0)
  # Malformed records stop with the same message as the R reader.
  writeLines("chr1\t100\t200\t1", files[1])
  expect_error(MscoreDMR:::.prepare_mscore_store(files, cpg, 1L, root), "six")
})

test_that("dense in-memory loci round-trip the long locus table", {
  dt <- data.table::data.table(sample = c("b", "a", "a", "b", "a"), chr = c(rep("chr1", 4), "chr2"),
    pos = c(200, 300, 100, 100, 50), S = c(1, 2, 3, 4, 5), T = c(5, 6, 7, 8, 9))
  dense <- MscoreDMR:::.dense_loci_from_long(dt, c("a", "b"))
  expect_equal(dense$chr1$pos, c(100, 200, 300))
  expect_equal(unname(dense$chr1$T), matrix(c(7, 0, 6, 8, 5, 0), 3))
  table <- MscoreDMR:::.store_locus_table(list(loci = dense))
  expected <- data.table::copy(dt)
  data.table::setorder(expected, chr, pos, sample)
  data.table::setcolorder(expected, c("sample", "chr", "pos", "T", "S"))
  expected[, M := S / T]
  expect_equal(table, expected)
  sums <- .Call(MscoreDMR:::C_region_sums, dense$chr1$pos, dense$chr1$S, dense$chr1$T,
                c(100, 101, 1), c(300, 300, 99))
  expect_equal(unname(sums$T), matrix(c(13, 6, 0, 13, 5, 0), 3))
})
