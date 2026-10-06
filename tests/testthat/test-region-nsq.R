brute_force_nsq <- function(reads, regions) {
  N <- numeric(nrow(regions))
  empty <- 0
  for (r in seq_len(nrow(regions))) {
    overlap <- reads$lo <= regions$end[r] & reads$hi >= regions$start[r]
    shared <- pmax(0, pmin(reads$b, regions$b[r]) - pmax(reads$a, regions$a[r]) + 1)
    N[r] <- sum(reads$count[overlap] * shared[overlap]^2)
    empty <- empty + sum(reads$nd[overlap & shared == 0])
  }
  list(N = N, empty = empty)
}

test_that("C region N_sq kernel matches a brute-force reference", {
  set.seed(42)
  positions <- sort(sample(1:5000, 600))
  for (iteration in 1:20) {
    n_reads <- sample(c(0, 1, 50, 400), 1)
    first <- sample(length(positions), n_reads, replace = TRUE)
    last <- pmin(length(positions), first + sample(0:25, n_reads, replace = TRUE))
    # Reads start at a CpG or slightly before it; ends extend past the last CpG.
    reads <- data.frame(lo = positions[first] - sample(0:3, n_reads, replace = TRUE),
                        hi = positions[last] + sample(0:3, n_reads, replace = TRUE))
    reads$a <- findInterval(reads$lo - 1, positions) + 1L
    reads$b <- findInterval(reads$hi, positions)
    reads$count <- as.numeric(sample(1:20, n_reads, replace = TRUE))
    reads$nd <- as.numeric(sample(1:3, n_reads, replace = TRUE))
    reads <- reads[order(reads$lo, reads$hi), ]
    start <- sample(1:5000, 30)
    end <- start + sample(0:800, 30, replace = TRUE)
    regions <- data.frame(start = as.numeric(start), end = as.numeric(end),
      a = findInterval(start - 1, positions) + 1L, b = findInterval(end, positions))
    observed <- .Call(MscoreDMR:::C_region_nsq, as.numeric(reads$lo),
      as.numeric(reads$hi), reads$a, reads$b, reads$count, reads$nd,
      regions$start, regions$end, regions$a, regions$b)
    expect_equal(observed, brute_force_nsq(reads, regions), tolerance = 0)
  }
})

test_that("C region N_sq kernel validates its inputs", {
  call <- function(lo, a = 1L) .Call(MscoreDMR:::C_region_nsq, lo, lo, a, a,
    1, 1, 1, 1, 1L, 1L)
  expect_error(call(1L), "double")
  expect_error(call(1, 1), "integer")
  expect_error(.Call(MscoreDMR:::C_region_nsq, c(2, 1), c(2, 1), 1:2, 1:2,
    c(1, 1), c(1, 1), 1, 1, 1L, 1L), "sorted")
})
