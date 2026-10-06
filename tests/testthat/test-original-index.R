test_that("original indexes preserve duplicates, crossing reads and exact matrices", {
  root <- tempfile("original-index-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  positions <- c(99, 100, 101, 199, 200, 201, 301)
  data.table::fwrite(data.frame(chr = "chr1", pos = positions), cpg,
                     sep = "\t", col.names = FALSE)
  # Identical rows are real contributions; duplicate retrieval is not.
  records <- data.frame(chr = "chr1", start = c(99,99,100,199,200,301,100),
    end = c(201,201,200,201,201,301,201),
    hap = c("100000","100000","0000","111","00","1","0X000"),
    count = c(2,3,5,7,11,2^32,4), strand = c("+","+","-","+","-","+","+"))
  records <- records[order(records$start, records$end), ]
  files <- setNames(file.path(root, paste0(c("a", "b"), ".tsv")), c("a", "b"))
  for (path in files) data.table::fwrite(records, path, sep = "\t", col.names = FALSE)
  legacy <- MscoreDMR:::.prepare_mscore_store(files, cpg, chunk_nrows = 2L,
    read_backend = "tabix")
  on.exit(unlink(legacy$root, recursive = TRUE), add = TRUE)
  packed <- vapply(files, Rsamtools::bgzip, character(1))
  for (i in seq_along(packed)) Rsamtools::indexTabix(packed[i], seq = 1,
    start = 2, end = 3, zeroBased = i == 2L)
  before <- tools::md5sum(c(packed, paste0(packed, ".tbi")))
  indexed <- MscoreDMR:::.prepare_mscore_store(packed, cpg, chunk_nrows = 2L,
    read_backend = "tabix")
  on.exit(unlink(indexed$root, recursive = TRUE), add = TRUE)
  expect_length(indexed$read_sources, 2L)
  expect_length(indexed$read_chunks, 0L)
  expect_length(list.files(indexed$root, pattern = "reads-"), 0L)
  expect_equal(indexed$qc, legacy$qc)
  expect_identical(indexed$loci, legacy$loci)
  regions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(
    c(99, 100, 101, 150, 199, 200, 201, 250, 301),
    c(301, 200, 201, 180, 199, 200, 201, 300, 301)))
  reference <- MscoreDMR:::.region_matrices_store(legacy, regions, names(files))
  observed <- MscoreDMR:::.region_matrices_store(indexed, regions, names(files))
  expect_equal(observed, reference, tolerance = 0)
  # The compact C backend must reproduce the tabix paths exactly, for plain
  # and indexed inputs alike (including the empty-hit QC count).
  for (input in list(files, packed)) {
    compact <- MscoreDMR:::.prepare_mscore_store(input, cpg, chunk_nrows = 2L)
    on.exit(unlink(compact$root, recursive = TRUE), add = TRUE)
    expect_identical(compact$read_backend, "compact")
    expect_length(compact$read_sources, 0L)
    expect_length(compact$read_chunks, 0L)
    expect_equal(compact$qc, legacy$qc)
    expect_equal(MscoreDMR:::.region_matrices_store(compact, regions, names(files)),
                 reference, tolerance = 0)
  }
  rt <- data.table::data.table(chr = "chr1", start = 99, end = 301, id = 1L)
  ann <- indexed$annotations[[1]]
  for (source in indexed$read_sources) {
    narrow <- MscoreDMR:::.original_region_contributions(source, rt, ann, 1L, window_bp = 100L)
    wide <- MscoreDMR:::.original_region_contributions(source, rt, ann, 100L, window_bp = 1000L)
    expect_equal(narrow, wide, tolerance = 0)
    expect_equal(narrow$N, 5 * 6^2 + 5 * 4^2 + 7 * 3^2 + 11 * 2^2 + 2^32)
  }
  expect_equal(MscoreDMR:::.region_matrices_store(indexed,
    GenomicRanges::GRanges(), names(files))$N_sq_matrix,
    matrix(numeric(), 0, 2, dimnames = list(NULL, names(files))))
  expect_identical(tools::md5sum(names(before)), before)
  mixed_files <- c(a = unname(packed[1]), b = unname(files[2]))
  mixed <- MscoreDMR:::.prepare_mscore_store(mixed_files, cpg, chunk_nrows = 2L,
    read_backend = "tabix")
  on.exit(unlink(mixed$root, recursive = TRUE), add = TRUE)
  expect_length(mixed$read_sources, 1L)
  expect_true(all(vapply(mixed$read_chunks, function(x) x$sample == "b", logical(1))))
  expect_equal(MscoreDMR:::.region_matrices_store(mixed, regions, names(files)),
               reference, tolerance = 0)
  # Detect changed inputs without rebuilding or editing the original index.
  changed <- indexed$read_sources[[1]]
  changed$fingerprint$size[1] <- changed$fingerprint$size[1] + 1
  expect_error(MscoreDMR:::.original_region_contributions(changed, rt, ann, 1L),
    "changed during analysis")
})

test_that("window ownership is per read-region pair, not per read", {
  root <- tempfile("window-index-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  cpg <- file.path(root, "cpg.tsv")
  data.table::fwrite(data.frame(chr = "chr1", pos = c(50, 100, 150, 200, 250)),
    cpg, sep = "\t", col.names = FALSE)
  file <- file.path(root, "reads.tsv")
  writeLines(c("chr1\t50\t250\t10000\t3\t+", "chr1\t50\t250\t10000\t4\t+"), file)
  packed <- Rsamtools::bgzip(file)
  Rsamtools::indexTabix(packed, seq = 1L, start = 2L, end = 3L, zeroBased = FALSE)
  source <- MscoreDMR:::.inspect_mhap_index(packed)
  source$sample <- "a"
  source$key <- "annotation1"
  rt <- data.table::data.table(chr = "chr1", start = c(50, 100, 151),
    end = c(250, 200, 251), id = 1:3)
  mapper <- getFromNamespace(".map_reads_to_cpg", "MscoreDMR")
  query_flags <- logical()
  testthat::local_mocked_bindings(.map_reads_to_cpg = function(raw_reads, cpg_cache,
                                                               compute_global = TRUE) {
    query_flags <<- c(query_flags, compute_global)
    mapper(raw_reads, cpg_cache, compute_global)
  }, .package = "MscoreDMR")
  result <- MscoreDMR:::.original_region_contributions(source, rt,
    list(chr1 = c(50, 100, 150, 200, 250)), chunk_nrows = 1L, window_bp = 100L)
  expect_named(result, c("N", "empty_hits"))
  expect_gt(length(query_flags), 0L)
  expect_false(any(query_flags))
  expect_equal(result$N, 7 * c(5, 3, 2)^2)
})

test_that("only N_sq queries skip state; initial single-site S/T remain exact", {
  raw <- data.table::data.table(sample = "a", annotation_key = "annotation1",
    chr = "chr1", read_start = c(100, 100, 100, 100), read_end = 300,
    original_hap_string = c("100", "000", "0X0", "00"),
    count = c(3, 2, 50, 60), strand = "+")
  cache <- list(annotation1 = list(chr1 = c(100, 200, 300)))
  initial <- MscoreDMR:::.map_reads_to_cpg(raw, cache)
  query <- MscoreDMR:::.map_reads_to_cpg(raw, cache, compute_global = FALSE)
  expect_true("is_methylated_global" %in% names(initial$reads))
  expect_false("is_methylated_global" %in% names(query$reads))
  expect_equal(initial$qc, query$qc)
  expect_equal(initial$qc$valid_mhap_records, 2L)
  expect_equal(initial$qc$invalid_mhap_records, 2L)
  expect_equal(query$reads, initial$reads[, !"is_methylated_global"])
  delta <- MscoreDMR:::.add_cpg_events(list(S = numeric(4), T = numeric(4)), initial$reads)
  expect_equal(head(cumsum(delta$S), -1L), rep(3, 3))
  expect_equal(head(cumsum(delta$T), -1L), rep(5, 3))
  # Even the locally unmethylated CpGs at 200 and 300 retain the global signal.
  empty <- MscoreDMR:::.map_reads_to_cpg(raw[0], cache, compute_global = FALSE)
  expect_equal(nrow(empty$reads), 0L)
  expect_false("is_methylated_global" %in% names(empty$reads))
})

test_that("broken or wrong-column indexes fail rather than silently rebuilding", {
  root <- tempfile("bad-index-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  file <- file.path(root, "x.tsv")
  writeLines("chr1\t100\t200\t10\t1\t+", file)
  packed <- Rsamtools::bgzip(file)
  writeLines("not an index", paste0(packed, ".tbi"))
  expect_error(MscoreDMR:::.inspect_mhap_index(packed), "Cannot use mHap index")
  unlink(paste0(packed, ".tbi"))
  Rsamtools::indexTabix(packed, seq = 1L, start = 2L, end = 2L)
  expect_error(MscoreDMR:::.inspect_mhap_index(packed), "columns 1/2/3")
  expect_null(MscoreDMR:::.inspect_mhap_index(file))
})

test_that("indexed samples use their own annotations across chromosomes", {
  root <- tempfile("multi-annotation-index-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  samples <- c("a", "b")
  files <- setNames(file.path(root, paste0(samples, ".tsv")), samples)
  annotations <- setNames(file.path(root, paste0(samples, ".cpg")), samples)
  for (i in seq_along(samples)) {
    pos <- if (i == 1L) c(10, 20, 30) else c(10, 25, 30)
    data.table::fwrite(data.frame(chr = rep(c("chr1", "chr2"), each = 3),
      pos = rep(pos, 2)), annotations[i], sep = "\t", col.names = FALSE)
    data.table::fwrite(data.frame(chr = c("chr1", "chr2"), start = 10,
      end = 30, string = c("100", "000"), count = c(2, 3), strand = "+"),
      files[i], sep = "\t", col.names = FALSE)
  }
  plain <- MscoreDMR:::.prepare_mscore_store(files, annotations, chunk_nrows = 1L,
    read_backend = "tabix")
  on.exit(unlink(plain$root, recursive = TRUE), add = TRUE)
  packed <- vapply(files, Rsamtools::bgzip, character(1))
  for (p in packed) Rsamtools::indexTabix(p, seq = 1L, start = 2L, end = 3L)
  indexed <- MscoreDMR:::.prepare_mscore_store(packed, annotations, chunk_nrows = 1L,
    read_backend = "tabix")
  on.exit(unlink(indexed$root, recursive = TRUE), add = TRUE)
  # Interleave chromosomes and include a region absent from the read index.
  regions <- GenomicRanges::GRanges(c("chr2", "chr1", "chr2", "chrAbsent"),
    IRanges::IRanges(c(20, 10, 10, 1), c(25, 30, 30, 100)))
  expect_equal(MscoreDMR:::.region_matrices_store(indexed, regions, samples),
    MscoreDMR:::.region_matrices_store(plain, regions, samples), tolerance = 0)
  expect_length(indexed$read_sources, 2L)
  expect_length(indexed$read_chunks, 0L)
  compact <- MscoreDMR:::.prepare_mscore_store(packed, annotations, chunk_nrows = 1L)
  on.exit(unlink(compact$root, recursive = TRUE), add = TRUE)
  expect_equal(MscoreDMR:::.region_matrices_store(compact, regions, samples),
    MscoreDMR:::.region_matrices_store(plain, regions, samples), tolerance = 0)
})
