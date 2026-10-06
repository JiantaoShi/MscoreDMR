# Internal data-ingestion and CpG-alignment helpers -------------------------

.as_path_vector <- function(paths, argument) {
  if (is.list(paths) && all(vapply(paths, function(path) {
    is.character(path) && length(path) == 1L
  }, logical(1)))) {
    paths <- unlist(paths, use.names = TRUE)
  }
  if (!is.character(paths) || !length(paths) || anyNA(paths) ||
      any(!nzchar(paths))) {
    stop("`", argument, "` must contain non-empty file paths.",
         call. = FALSE)
  }
  paths
}

.normalise_mhap_files <- function(mhap_files) {
  mhap_files <- .as_path_vector(mhap_files, "mhap_files")
  if (length(mhap_files) < 2L) {
    stop("`mhap_files` must contain at least two samples.", call. = FALSE)
  }
  if (any(!file.exists(mhap_files))) {
    stop(
      "mHap file(s) not found: ",
      paste(mhap_files[!file.exists(mhap_files)], collapse = ", "),
      call. = FALSE
    )
  }
  sample_names <- names(mhap_files)
  if (is.null(sample_names) || any(!nzchar(sample_names))) {
    sample_names <- paste0("sample", seq_along(mhap_files))
  }
  if (anyDuplicated(sample_names)) {
    stop("Sample names derived from `mhap_files` must be unique.",
         call. = FALSE)
  }
  names(mhap_files) <- sample_names
  mhap_files
}

.normalise_cpg_files <- function(cpg_files, sample_names) {
  cpg_files <- .as_path_vector(cpg_files, "cpg_files")
  if (length(cpg_files) == 1L) {
    cpg_files <- rep(cpg_files, length(sample_names))
    names(cpg_files) <- sample_names
  } else {
    if (length(cpg_files) != length(sample_names)) {
      stop("`cpg_files` must contain one shared path or one path per sample.",
           call. = FALSE)
    }
    cpg_names <- names(cpg_files)
    if (is.null(cpg_names) || any(!nzchar(cpg_names)) ||
        anyDuplicated(cpg_names)) {
      stop("Multiple `cpg_files` must have unique, non-empty sample names.",
           call. = FALSE)
    }
    if (!setequal(cpg_names, sample_names)) {
      stop("Names of `cpg_files` must match names of `mhap_files`.",
           call. = FALSE)
    }
    cpg_files <- cpg_files[match(sample_names, cpg_names)]
  }
  if (any(!file.exists(cpg_files))) {
    stop(
      "CpG annotation file(s) not found: ",
      paste(unique(cpg_files[!file.exists(cpg_files)]), collapse = ", "),
      call. = FALSE
    )
  }
  normalised <- normalizePath(cpg_files, winslash = "/", mustWork = TRUE)
  names(normalised) <- sample_names
  normalised
}

.read_mhap_file <- function(path) {
  data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    select = 1:6,
    col.names = c("chr", "read_start", "read_end", "original_hap_string",
                  "count", "strand"),
    colClasses = c(
      "character", "numeric", "numeric", "character", "numeric", "character"
    ),
    showProgress = FALSE
  )
}

.read_cpg_file <- function(path, spans) {
  if (file.exists(paste0(path, ".tbi"))) {
    pieces <- vector("list", nrow(spans))
    for (index in seq_len(nrow(spans))) {
      query <- GenomicRanges::GRanges(
        spans$chr[index],
        IRanges::IRanges(spans$span_start[index], spans$span_end[index])
      )
      lines <- Rsamtools::scanTabix(path, param = query)[[1L]]
      if (!length(lines)) {
        pieces[[index]] <- data.table::data.table(
          chr = character(), pos = numeric()
        )
      } else {
        fields <- strsplit(lines, "\t", fixed = TRUE)
        pieces[[index]] <- data.table::data.table(
          chr = vapply(fields, `[`, character(1L), 1L),
          pos = suppressWarnings(as.numeric(vapply(
            fields, `[`, character(1L), 2L
          )))
        )
      }
    }
    return(data.table::rbindlist(pieces, use.names = TRUE))
  }

  annotation <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    select = 1:2,
    col.names = c("chr", "pos"),
    colClasses = list(character = 1L, numeric = 2L),
    showProgress = FALSE
  )
  if (!nrow(annotation) || !nrow(spans)) return(annotation[0])
  pieces <- lapply(seq_len(nrow(spans)), function(index) {
    annotation[
      chr == spans$chr[index] &
        pos >= spans$span_start[index] & pos <= spans$span_end[index]
    ]
  })
  data.table::rbindlist(pieces, use.names = TRUE)
}

.load_raw_mhap_reads <- function(mhap_files, annotation_by_sample) {
  pieces <- vector("list", length(mhap_files))
  for (sample_index in seq_along(mhap_files)) {
    sample_name <- names(mhap_files)[sample_index]
    message("Reading sample ", sample_name, "...")
    piece <- .read_mhap_file(mhap_files[[sample_index]])
    piece[, `:=`(
      sample = sample_name,
      annotation_key = unname(annotation_by_sample[sample_name])
    )]
    pieces[[sample_index]] <- piece
  }
  data.table::rbindlist(pieces, use.names = TRUE)
}

.load_cpg_cache <- function(annotation_paths, raw_reads) {
  cache <- stats::setNames(vector("list", length(annotation_paths)),
                           names(annotation_paths))
  for (key_value in names(annotation_paths)) {
    relevant <- raw_reads[
      annotation_key == key_value &
        is.finite(read_start) & is.finite(read_end) &
        read_start <= read_end & !is.na(chr) & nzchar(chr)
    ]
    spans <- relevant[, .(
      span_start = min(read_start),
      span_end = max(read_end)
    ), by = chr]
    message("Loading CpG annotation ", key_value, "...")
    annotation <- .read_cpg_file(annotation_paths[[key_value]], spans)
    annotation <- annotation[
      !is.na(chr) & nzchar(chr) & is.finite(pos) & pos >= 1
    ]
    annotation[, pos := as.numeric(pos)]
    annotation <- unique(annotation, by = c("chr", "pos"))
    data.table::setorder(annotation, chr, pos)
    cache[[key_value]] <- split(annotation$pos, annotation$chr)
  }
  cache
}

.map_reads_to_cpg <- function(raw_reads, cpg_cache, compute_global = TRUE) {
  total_records <- nrow(raw_reads)
  invalid_string <- is.na(raw_reads$original_hap_string) |
    !grepl("^[01]+$", raw_reads$original_hap_string)
  invalid_basic <- is.na(raw_reads$chr) | !nzchar(raw_reads$chr) |
    !is.finite(raw_reads$read_start) | !is.finite(raw_reads$read_end) |
    raw_reads$read_start > raw_reads$read_end |
    !is.finite(raw_reads$count) | raw_reads$count <= 0 |
    is.na(raw_reads$strand) | !raw_reads$strand %in% c("+", "-") |
    invalid_string

  valid_basic <- data.table::copy(raw_reads[!invalid_basic])
  if (nrow(valid_basic)) {
    valid_basic <- valid_basic[, .(
      count = sum(count),
      source_record_count = .N
    ), by = .(
      sample, annotation_key, chr, read_start, read_end,
      original_hap_string, strand
    )]
  }
  # Initial locus counting needs the full-read state. N_sq-only region queries
  # reuse validation/mapping but do not need to compute or store that state.
  if (compute_global) {
    valid_basic[, is_methylated_global := grepl(
      "1", original_hap_string, fixed = TRUE)]
  }
  valid_basic[, `:=`(
    cpg_index_start = NA_integer_,
    cpg_index_end = NA_integer_,
    mapping_valid = FALSE,
    missing_annotation = FALSE,
    length_mismatch = FALSE
  )]

  if (nrow(valid_basic)) {
    groups <- unique(valid_basic[, .(annotation_key, chr)])
    for (group_index in seq_len(nrow(groups))) {
      key <- groups$annotation_key[group_index]
      chromosome <- groups$chr[group_index]
      row_index <- which(
        valid_basic$annotation_key == key & valid_basic$chr == chromosome
      )
      positions <- cpg_cache[[key]][[chromosome]]
      if (is.null(positions) || !length(positions)) {
        valid_basic$missing_annotation[row_index] <- TRUE
        next
      }
      first_index <- findInterval(
        valid_basic$read_start[row_index] - 1, positions
      ) + 1L
      last_index <- findInterval(valid_basic$read_end[row_index], positions)
      covered <- pmax(0L, last_index - first_index + 1L)
      expected <- nchar(valid_basic$original_hap_string[row_index])
      valid_in_group <- covered > 0L & covered == expected
      valid_basic$cpg_index_start[row_index] <- first_index
      valid_basic$cpg_index_end[row_index] <- last_index
      valid_basic$missing_annotation[row_index] <- covered == 0L
      valid_basic$length_mismatch[row_index] <- covered > 0L & covered != expected
      valid_basic$mapping_valid[row_index] <- valid_in_group
    }
  }

  mismatch_records <- if (nrow(valid_basic)) {
    sum(valid_basic$source_record_count[valid_basic$length_mismatch])
  } else 0L
  missing_records <- if (nrow(valid_basic)) {
    sum(valid_basic$source_record_count[valid_basic$missing_annotation])
  } else 0L
  mapped_reads <- valid_basic[which(valid_basic$mapping_valid)]
  valid_record_count <- if (nrow(mapped_reads)) {
    sum(mapped_reads$source_record_count)
  } else 0L
  qc <- list(
    total_mhap_records = as.integer(total_records),
    valid_mhap_records = as.integer(valid_record_count),
    invalid_mhap_records = as.integer(total_records - valid_record_count),
    cpg_length_mismatch_records = as.integer(mismatch_records),
    invalid_haplotype_string_records = as.integer(sum(invalid_string)),
    missing_cpg_annotation_records = as.integer(missing_records),
    reads_without_local_cpg = 0L,
    regions_dropped_incomplete_coverage = 0L
  )
  if (qc$invalid_mhap_records > 0L) {
    message(
      "Discarded ", qc$invalid_mhap_records, " invalid mHap record(s): ",
      qc$cpg_length_mismatch_records, " CpG/string length mismatch, ",
      qc$invalid_haplotype_string_records, " invalid string, ",
      qc$missing_cpg_annotation_records, " without mapped CpGs."
    )
  }
  data.table::setorder(mapped_reads, annotation_key, chr, read_start, sample)
  list(reads = mapped_reads, qc = qc)
}

.build_master_loci <- function(master_read_dt, cpg_cache) {
  if (!nrow(master_read_dt)) {
    return(data.table::data.table(
      sample = character(), chr = character(), pos = numeric(),
      T = numeric(), S = numeric(), M = numeric()
    ))
  }

  max_expanded_rows <- 2000000L
  buffer_table_limit <- 8L
  buffer_row_limit <- max_expanded_rows
  read_groups <- unique(master_read_dt[, .(annotation_key, chr)])
  pieces <- vector("list", nrow(read_groups))

  for (group_index in seq_len(nrow(read_groups))) {
    key <- read_groups$annotation_key[group_index]
    chromosome <- read_groups$chr[group_index]
    reads <- master_read_dt[
      annotation_key == key & chr == chromosome,
      .(
        sample, cpg_index_start, cpg_index_end, count,
        is_methylated_global
      )
    ]
    positions <- cpg_cache[[key]][[chromosome]]

    message(
      "Building master loci for ", key, " / ", chromosome, ": ",
      format(nrow(reads), big.mark = ",", scientific = FALSE),
      " reads..."
    )

    if (!nrow(reads) || is.null(positions) || !length(positions)) {
      message("No mapped CpGs for ", key, " / ", chromosome, ".")
      next
    }

    ncpg <- as.double(
      reads$cpg_index_end - reads$cpg_index_start + 1L
    )
    valid_ncpg <- is.finite(ncpg) & ncpg > 0
    if (any(!valid_ncpg)) {
      message(
        "Skipping ", sum(!valid_ncpg),
        " read(s) with invalid CpG spans in ", key, " / ", chromosome,
        "."
      )
      reads <- reads[valid_ncpg]
      ncpg <- ncpg[valid_ncpg]
    }
    if (!nrow(reads)) {
      message("No valid CpG spans for ", key, " / ", chromosome, ".")
      next
    }

    cumulative_ncpg <- c(0, cumsum(ncpg))
    total_chunks <- 0
    scan_start <- 1L
    while (scan_start <= nrow(reads)) {
      if (ncpg[scan_start] > max_expanded_rows) {
        total_chunks <- total_chunks +
          ceiling(ncpg[scan_start] / max_expanded_rows)
        scan_start <- scan_start + 1L
      } else {
        cumulative_limit <- cumulative_ncpg[scan_start] +
          max_expanded_rows
        scan_end <- min(
          nrow(reads),
          findInterval(cumulative_limit, cumulative_ncpg) - 1L
        )
        scan_end <- max(scan_start, scan_end)
        total_chunks <- total_chunks + 1
        scan_start <- scan_end + 1L
      }
    }

    message(
      "Processing ", format(total_chunks, big.mark = ","),
      " chunk(s) for ", key, " / ", chromosome, "."
    )

    aggregation_buffer <- vector("list", 0L)
    buffered_rows <- 0L
    chromosome_accumulator <- NULL

    flush_buffer <- function(buffer, accumulator) {
      if (!length(buffer)) {
        return(accumulator)
      }
      combined <- data.table::rbindlist(buffer, use.names = TRUE)
      combined <- combined[, .(
        T = sum(as.numeric(T)),
        S = sum(as.numeric(S))
      ), by = .(sample, pos)]
      if (is.null(accumulator)) {
        return(combined)
      }
      data.table::rbindlist(
        list(accumulator, combined), use.names = TRUE
      )[, .(
        T = sum(as.numeric(T)),
        S = sum(as.numeric(S))
      ), by = .(sample, pos)]
    }

    add_to_buffer <- function(chunk_result, buffer, buffered_rows,
                              accumulator) {
      if (length(buffer) &&
          (length(buffer) >= buffer_table_limit ||
           buffered_rows + nrow(chunk_result) > buffer_row_limit)) {
        accumulator <- flush_buffer(buffer, accumulator)
        buffer <- vector("list", 0L)
        buffered_rows <- 0L
      }

      buffer[[length(buffer) + 1L]] <- chunk_result
      buffered_rows <- buffered_rows + nrow(chunk_result)

      if (length(buffer) >= buffer_table_limit ||
          buffered_rows >= buffer_row_limit) {
        accumulator <- flush_buffer(buffer, accumulator)
        buffer <- vector("list", 0L)
        buffered_rows <- 0L
      }

      list(
        buffer = buffer,
        buffered_rows = buffered_rows,
        accumulator = accumulator
      )
    }

    chunk_number <- 0
    read_start_index <- 1L
    while (read_start_index <= nrow(reads)) {
      if (ncpg[read_start_index] > max_expanded_rows) {
        remaining <- ncpg[read_start_index]
        offset <- 0
        while (remaining > 0) {
          subrange_length <- min(
            as.double(max_expanded_rows), remaining
          )
          subrange_start <-
            as.double(reads$cpg_index_start[read_start_index]) + offset
          subrange_end <- subrange_start + subrange_length - 1
          cpg_idx <- seq.int(subrange_start, subrange_end)
          count_numeric <- as.numeric(reads$count[read_start_index])
          add_S <- count_numeric *
            as.numeric(reads$is_methylated_global[read_start_index])
          expanded <- data.table::data.table(
            sample = reads$sample[read_start_index],
            pos = positions[cpg_idx],
            T = count_numeric,
            S = add_S
          )
          chunk_result <- expanded[, .(
            T = sum(as.numeric(T)),
            S = sum(as.numeric(S))
          ), by = .(sample, pos)]

          chunk_number <- chunk_number + 1
          buffered <- add_to_buffer(
            chunk_result, aggregation_buffer, buffered_rows,
            chromosome_accumulator
          )
          aggregation_buffer <- buffered$buffer
          buffered_rows <- buffered$buffered_rows
          chromosome_accumulator <- buffered$accumulator

          if (chunk_number == 1 || chunk_number == total_chunks ||
              chunk_number %% 10 == 0) {
            message(
              "Processing chunk ", chunk_number, " of ", total_chunks,
              " for ", key, " / ", chromosome, "..."
            )
          }

          remaining <- remaining - subrange_length
          offset <- offset + subrange_length
          rm(
            cpg_idx, count_numeric, add_S, expanded, chunk_result,
            buffered
          )
          gc(verbose = FALSE)
        }
        read_start_index <- read_start_index + 1L
        next
      }

      cumulative_limit <- cumulative_ncpg[read_start_index] +
        max_expanded_rows
      read_end_index <- min(
        nrow(reads),
        findInterval(cumulative_limit, cumulative_ncpg) - 1L
      )
      read_end_index <- max(read_start_index, read_end_index)
      chunk_rows <- seq.int(read_start_index, read_end_index)
      chunk_ncpg <- as.integer(ncpg[chunk_rows])
      read_row <- rep.int(seq_along(chunk_rows), times = chunk_ncpg)
      chunk_cpg_start <- reads$cpg_index_start[chunk_rows]
      cpg_idx <- chunk_cpg_start[read_row] + sequence(chunk_ncpg) - 1L
      chunk_sample <- reads$sample[chunk_rows]
      count_numeric <- as.numeric(reads$count[chunk_rows])
      add_S <- count_numeric *
        as.numeric(reads$is_methylated_global[chunk_rows])
      expanded <- data.table::data.table(
        sample = chunk_sample[read_row],
        pos = positions[cpg_idx],
        T = count_numeric[read_row],
        S = add_S[read_row]
      )
      chunk_result <- expanded[, .(
        T = sum(as.numeric(T)),
        S = sum(as.numeric(S))
      ), by = .(sample, pos)]

      chunk_number <- chunk_number + 1
      buffered <- add_to_buffer(
        chunk_result, aggregation_buffer, buffered_rows,
        chromosome_accumulator
      )
      aggregation_buffer <- buffered$buffer
      buffered_rows <- buffered$buffered_rows
      chromosome_accumulator <- buffered$accumulator

      if (chunk_number == 1 || chunk_number == total_chunks ||
          chunk_number %% 10 == 0) {
        message(
          "Processing chunk ", chunk_number, " of ", total_chunks,
          " for ", key, " / ", chromosome, "..."
        )
      }

      read_start_index <- read_end_index + 1L
      rm(
        chunk_rows, chunk_ncpg, read_row, chunk_cpg_start, cpg_idx,
        chunk_sample, count_numeric, add_S, expanded, chunk_result,
        buffered
      )
      gc(verbose = FALSE)
    }

    chromosome_accumulator <- flush_buffer(
      aggregation_buffer, chromosome_accumulator
    )
    if (!is.null(chromosome_accumulator) &&
        nrow(chromosome_accumulator)) {
      chromosome_accumulator[, chr := chromosome]
      data.table::setcolorder(
        chromosome_accumulator, c("sample", "chr", "pos", "T", "S")
      )
      pieces[[group_index]] <- chromosome_accumulator
    }
    message("Completed ", key, " / ", chromosome, ".")
    rm(
      reads, positions, ncpg, valid_ncpg, cumulative_ncpg,
      aggregation_buffer, chromosome_accumulator, flush_buffer,
      add_to_buffer
    )
    gc(verbose = FALSE)
  }

  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) {
    return(data.table::data.table(
      sample = character(), chr = character(), pos = numeric(),
      T = numeric(), S = numeric(), M = numeric()
    ))
  }
  loci <- data.table::rbindlist(pieces, use.names = TRUE)
  loci <- loci[, .(
    T = sum(as.numeric(T)),
    S = sum(as.numeric(S))
  ), by = .(sample, chr, pos)]
  loci[, M := S / T]
  data.table::setorder(loci, chr, pos, sample)
  loci
}

.prepare_mscore_data <- function(mhap_files, cpg_files) {
  mhap_files <- .normalise_mhap_files(mhap_files)
  sample_names <- names(mhap_files)
  cpg_files <- .normalise_cpg_files(cpg_files, sample_names)
  unique_paths <- unique(unname(cpg_files))
  annotation_paths <- stats::setNames(
    unique_paths, paste0("annotation", seq_along(unique_paths))
  )
  path_to_key <- stats::setNames(names(annotation_paths), annotation_paths)
  annotation_by_sample <- stats::setNames(
    unname(path_to_key[cpg_files]), sample_names
  )

  raw_reads <- .load_raw_mhap_reads(mhap_files, annotation_by_sample)
  cpg_cache <- .load_cpg_cache(annotation_paths, raw_reads)
  mapped <- .map_reads_to_cpg(raw_reads, cpg_cache)
  master_locus_dt <- .build_master_loci(mapped$reads, cpg_cache)
  list(
    mhap_files = mhap_files,
    cpg_files = cpg_files,
    annotation_paths = annotation_paths,
    annotation_by_sample = annotation_by_sample,
    cpg_cache = cpg_cache,
    master_read_dt = mapped$reads,
    master_locus_dt = master_locus_dt,
    qc = mapped$qc
  )
}

# Streaming production path. The helpers above remain available for small-data
# diagnostics and regression comparisons; exported analyses use this store.
.validate_chunk_nrows <- function(chunk_nrows) {
  if (length(chunk_nrows) != 1L || !is.numeric(chunk_nrows) ||
      !is.finite(chunk_nrows) || chunk_nrows < 1 ||
      chunk_nrows > .Machine$integer.max || chunk_nrows != floor(chunk_nrows)) {
    stop("`chunk_nrows` must be a positive integer.", call. = FALSE)
  }
  as.integer(chunk_nrows)
}

.stream_text_batches <- function(path, chunk_nrows, FUN) {
  # gzfile also reads uncompressed text. Keep the connection open: never use
  # fread(skip=...) repeatedly on compressed files.
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  repeat {
    lines <- readLines(connection, n = chunk_nrows, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines)) FUN(lines)
  }
  invisible(NULL)
}

.parse_mhap_lines <- function(lines) {
  if (any(lengths(strsplit(lines, "\t", fixed = TRUE)) < 6L)) {
    stop("mHap records must have at least six tab-separated fields.")
  }
  data.table::fread(text = paste0(paste(lines, collapse = "\n"), "\n"),
    sep = "\t", header = FALSE, select = 1:6,
    col.names = c("chr", "read_start", "read_end", "original_hap_string",
                  "count", "strand"),
    colClasses = c("character", "numeric", "numeric", "character",
                   "numeric", "character"), showProgress = FALSE)
}

.stream_cpg_annotation <- function(path, chunk_nrows) {
  pieces <- list()
  .stream_text_batches(path, chunk_nrows, function(lines) {
    dt <- data.table::fread(
      text = paste0(paste(lines, collapse = "\n"), "\n"),
      sep = "\t", header = FALSE, select = 1:2,
      col.names = c("chr", "pos"),
      colClasses = list(character = 1L, numeric = 2L), showProgress = FALSE)
    dt <- dt[!is.na(chr) & nzchar(chr) & is.finite(pos) & pos >= 1]
    pieces[[length(pieces) + 1L]] <<- dt
  })
  if (!length(pieces)) return(list())
  dt <- unique(data.table::rbindlist(pieces), by = c("chr", "pos"))
  data.table::setorder(dt, chr, pos)
  split(as.numeric(dt$pos), dt$chr)
}

.add_cpg_events <- function(delta, reads) {
  if (!nrow(reads)) return(delta)
  n <- nrow(reads)
  events <- data.table::data.table(
    index = c(reads$cpg_index_start, reads$cpg_index_end + 1L),
    T = rep(as.numeric(reads$count), 2L) * rep(c(1, -1), each = n),
    S = rep(as.numeric(reads$count) * as.numeric(reads$is_methylated_global),
            2L) * rep(c(1, -1), each = n))
  events <- events[, .(T = sum(T), S = sum(S)), by = index]
  # indices are unique HERE, including cancellations between starts and ends.
  delta$T[events$index] <- delta$T[events$index] + events$T
  delta$S[events$index] <- delta$S[events$index] + events$S
  delta
}

.inspect_mhap_index <- function(path) {
  index <- paste0(path, ".tbi")
  if (!file.exists(index)) return(NULL)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  index <- normalizePath(index, winslash = "/", mustWork = TRUE)
  tryCatch({
    tab <- Rsamtools::TabixFile(path, index = index)
    open(tab)
    on.exit(close(tab), add = TRUE)
    header <- Rsamtools::headerTabix(tab)
    if (!identical(as.integer(header$indexColumns), c(1L, 2L, 3L)) ||
        header$skip != 0L) {
      stop("Expected a headerless mHap index with chromosome/start/end columns 1/2/3.")
    }
    list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
      index = normalizePath(index, winslash = "/", mustWork = TRUE),
      chromosomes = header$seqnames,
      fingerprint = file.info(c(path, index))[, c("size", "mtime"), drop = FALSE])
  }, error = function(e) stop("Cannot use mHap index for ", path, ": ",
    conditionMessage(e), " Fix or explicitly remove the unusable index; it was not modified.",
    call. = FALSE))
}

.prepare_mscore_store <- function(mhap_files, cpg_files,
                                  chunk_nrows = 100000L, cache_dir = NULL,
                                  keep_reads = TRUE,
                                  read_backend = c("compact", "tabix")) {
  # "compact": validated reads are cached once as sorted per-chromosome arrays
  # and region N_sq is computed in C. "tabix": the original R implementation
  # that re-queries indexed reads for every region set (kept for validation).
  read_backend <- match.arg(read_backend)
  compact <- keep_reads && read_backend == "compact"
  # The C reader serves every path except the tabix validation backend.
  use_c_reader <- !(keep_reads && read_backend == "tabix")
  chunk_nrows <- .validate_chunk_nrows(chunk_nrows)
  mhap_files <- .normalise_mhap_files(mhap_files)
  samples <- names(mhap_files)
  cpg_files <- .normalise_cpg_files(cpg_files, samples)
  if (is.null(cache_dir)) cache_dir <- tempdir()
  if (!dir.exists(cache_dir) && !dir.create(cache_dir, recursive = TRUE)) {
    stop("Cannot create cache directory: ", cache_dir)
  }
  # A fresh private directory every time; never reuse stale caches or overwrite
  # user files. tempfile does not consume the permutation RNG stream.
  root <- tempfile("mscore-store-", tmpdir = cache_dir)
  dir.create(root)
  complete <- FALSE
  on.exit(if (!complete) unlink(root, recursive = TRUE), add = TRUE)
  paths <- unique(unname(cpg_files))
  annotation_paths <- stats::setNames(paths, paste0("annotation", seq_along(paths)))
  keys <- stats::setNames(names(annotation_paths)[match(cpg_files, paths)], samples)
  qc_names <- c("total_mhap_records", "valid_mhap_records", "invalid_mhap_records",
    "cpg_length_mismatch_records", "invalid_haplotype_string_records",
    "missing_cpg_annotation_records", "reads_without_local_cpg",
    "regions_dropped_incomplete_coverage")
  qc <- as.list(stats::setNames(rep(0, length(qc_names)), qc_names))
  # Covered loci per chromosome and sample: list(pos, T, S) with T > 0.
  sample_loci <- list()
  read_chunks <- list()
  read_sources <- list()
  read_paths <- list()
  annotations <- list()
  chunk_id <- 0L
  # Everything stays in memory except the binary read cache (and, for the
  # tabix validation backend, its indexed read shards).
  for (key in names(annotation_paths)) {
    message("Loading CpG coordinate index ", key, "...")
    positions <- .stream_cpg_annotation(annotation_paths[[key]], chunk_nrows)
    annotations[[key]] <- positions
    for (sample_name in samples[keys == key]) {
      source <- if (keep_reads && !compact) {
        .inspect_mhap_index(mhap_files[[sample_name]])
      } else NULL
      if (!is.null(source)) {
        source$sample <- sample_name
        source$key <- key
        read_sources[[sample_name]] <- source
        message("Reusing original mHap index for ", sample_name, ".")
      }
      deltas <- new.env(parent = emptyenv())
      batches <- 0L
      read_pieces <- list()
      scanned <- if (use_c_reader) {
        message("Streaming sample ", sample_name, " (C reader)...")
        .scan_mhap_file(mhap_files[[sample_name]], positions,
          if (compact) file.path(root, paste0("reads-", match(sample_name, samples))))
      }
      if (!is.null(scanned)) {
        for (nm in names(scanned$qc)) qc[[nm]] <- qc[[nm]] + scanned$qc[[nm]]
        for (chromosome in names(scanned$deltas)) {
          assign(chromosome, scanned$deltas[[chromosome]], deltas)
        }
        if (compact) read_paths[[sample_name]] <- as.list(scanned$read_files)
        rm(scanned)
      } else {
      message("Streaming sample ", sample_name, " (", chunk_nrows, " records/batch)...")
      .stream_text_batches(mhap_files[[sample_name]], chunk_nrows, function(lines) {
        batches <<- batches + 1L
        raw <- .parse_mhap_lines(lines)
        raw[, `:=`(sample = sample_name, annotation_key = key)]
        mapped <- .map_reads_to_cpg(raw, stats::setNames(list(positions), key))
        for (nm in qc_names) qc[[nm]] <<- qc[[nm]] + mapped$qc[[nm]]
        reads <- mapped$reads
        if (!is.null(source) && any(!unique(reads$chr) %in% source$chromosomes)) {
          stop("The original index omits a chromosome with valid reads: ", sample_name)
        }
        for (chromosome in unique(reads$chr)) {
          part <- reads[chr == chromosome]
          if (!exists(chromosome, deltas, inherits = FALSE)) {
            assign(chromosome, list(T = numeric(length(positions[[chromosome]]) + 1L),
              S = numeric(length(positions[[chromosome]]) + 1L)), deltas)
          }
          assign(chromosome, .add_cpg_events(get(chromosome, deltas), part), deltas)
        }
        if (compact && nrow(reads)) {
          for (chromosome in unique(reads$chr)) {
            piece <- reads[chr == chromosome, .(lo = read_start, hi = read_end,
              a = cpg_index_start, b = cpg_index_end, hap = original_hap_string,
              strand, count)]
            read_pieces[[chromosome]] <<- c(read_pieces[[chromosome]], list(piece))
          }
        }
        if (keep_reads && !compact && is.null(source) && nrow(reads)) {
          # Each indexed file is <= chunk_nrows validated records. This bounds
          # a single query even for very wide regions and supports unsorted TSV
          # inputs without sorting a whole sample in memory.
          chunk_id <<- chunk_id + 1L
          out <- reads[, .(chr, read_start, read_end, cpg_index_start,
            cpg_index_end, count, is_methylated_global,
            original_hap_string, strand)]
          data.table::setorder(out, chr, read_start, read_end)
          plain <- file.path(root, paste0("reads-", chunk_id, ".tsv"))
          # Tabix coordinate columns must be decimal integers, not 1e+08.
          data.table::fwrite(out, plain, sep = "\t", col.names = FALSE, scipen = 999)
          packed <- Rsamtools::bgzip(plain, overwrite = FALSE)
          Rsamtools::indexTabix(packed, seq = 1L, start = 2L, end = 3L,
                               zeroBased = FALSE)
          unlink(plain)
          spans <- out[, .(lo = min(read_start), hi = max(read_end)), by = chr]
          read_chunks[[length(read_chunks) + 1L]] <<- list(
            path = packed, sample = sample_name, key = key, spans = spans)
        }
        if (batches == 1L || batches %% 25L == 0L) {
          message("  ", sample_name, ": batch ", batches)
        }
      })
      if (compact) {
        read_paths[[sample_name]] <- .consolidate_read_pieces(
          read_pieces, root, match(sample_name, samples))
      }
      }
      for (chromosome in ls(deltas, all.names = TRUE)) {
        delta <- get(chromosome, deltas)
        T <- utils::head(cumsum(delta$T), -1L)
        S <- utils::head(cumsum(delta$S), -1L)
        if (any(!is.finite(T) | !is.finite(S) | T < 0 | S < 0 | S > T)) {
          stop("Invalid difference-array counts for ", sample_name, " / ", chromosome)
        }
        covered <- which(T > 0)
        sample_loci[[chromosome]][[sample_name]] <- list(
          pos = positions[[chromosome]][covered], T = T[covered], S = S[covered])
        rm(list = chromosome, envir = deltas)
      }
      rm(deltas)
      gc(verbose = FALSE)
    }
    rm(positions)
    gc(verbose = FALSE)
  }
  # Dense (locus x sample) S/T matrices per chromosome; zero = not covered.
  loci <- list()
  for (chromosome in sort(names(sample_loci))) {
    pieces <- sample_loci[[chromosome]]
    pos <- sort(unique(unlist(lapply(pieces, `[[`, "pos"), use.names = FALSE)))
    S <- matrix(0, length(pos), length(samples), dimnames = list(NULL, samples))
    T <- S
    for (sample_name in names(pieces)) {
      rows <- match(pieces[[sample_name]]$pos, pos)
      S[rows, sample_name] <- pieces[[sample_name]]$S
      T[rows, sample_name] <- pieces[[sample_name]]$T
    }
    loci[[chromosome]] <- list(pos = pos, S = S, T = T)
    sample_loci[[chromosome]] <- NULL
  }
  store <- list(root = root, mhap_files = mhap_files, cpg_files = cpg_files,
    sample_names = samples, annotations = annotations, annotation_by_sample = keys,
    loci = loci, read_chunks = read_chunks, read_sources = read_sources,
    read_backend = if (keep_reads) read_backend else "none",
    read_paths = read_paths, qc = qc, chunk_nrows = chunk_nrows)
  complete <- TRUE
  store
}

.consolidate_read_pieces <- function(read_pieces, root, sample_index) {
  # Merge one sample's in-memory batch pieces per chromosome. Identical mHap records
  # (also across batches) are merged first so that `nd` counts distinct
  # records, then rows sharing the same span are collapsed for N_sq.
  paths <- stats::setNames(character(length(read_pieces)), names(read_pieces))
  for (i in seq_along(read_pieces)) {
    dt <- data.table::rbindlist(read_pieces[[i]])
    dt <- dt[, .(count = sum(count)), by = .(lo, hi, a, b, hap, strand)]
    dt <- dt[, .(count = sum(count), nd = as.numeric(.N)), by = .(lo, hi, a, b)]
    data.table::setorder(dt, lo, hi)
    paths[i] <- file.path(root, paste0("reads-r-", sample_index, "-", i, ".bin"))
    .write_compact_reads(paths[i], dt)
  }
  as.list(paths)
}

# Binary layout shared with the C reader: n, lo, hi, a, b, count, nd.
.write_compact_reads <- function(path, dt) {
  connection <- file(path, "wb")
  on.exit(close(connection))
  writeBin(as.numeric(nrow(dt)), connection)
  writeBin(as.numeric(dt$lo), connection)
  writeBin(as.numeric(dt$hi), connection)
  writeBin(as.integer(dt$a), connection, size = 4L)
  writeBin(as.integer(dt$b), connection, size = 4L)
  writeBin(as.numeric(dt$count), connection)
  writeBin(as.numeric(dt$nd), connection)
  invisible(path)
}

.read_compact_reads <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  n <- readBin(connection, "double", 1L)
  list(lo = readBin(connection, "double", n), hi = readBin(connection, "double", n),
    a = readBin(connection, "integer", n, size = 4L),
    b = readBin(connection, "integer", n, size = 4L),
    count = readBin(connection, "double", n), nd = readBin(connection, "double", n))
}

# Run the C reader on one sample. Returns NULL when the read cache was
# requested but the input is not sorted by chromosome blocks and start; the
# caller then uses the R streaming reader for that sample.
.scan_mhap_file <- function(path, positions, read_prefix = NULL) {
  result <- .Call(C_scan_mhap, normalizePath(path, winslash = "/", mustWork = TRUE),
    as.character(names(positions)), unname(lapply(positions, as.numeric)),
    if (is.null(read_prefix)) NULL else as.character(read_prefix))
  if (!is.null(result$error)) {
    unlink(result$read_files)
    stop(result$error, " File: ", path, call. = FALSE)
  }
  if (result$unsorted) {
    unlink(result$read_files)
    message("  ", basename(path), " is not sorted by chromosome and start; ",
            "using the R reader for this sample.")
    return(NULL)
  }
  qc <- result$qc
  qc <- list(total_mhap_records = qc[1L], valid_mhap_records = qc[2L],
    invalid_mhap_records = qc[1L] - qc[2L], cpg_length_mismatch_records = qc[4L],
    invalid_haplotype_string_records = qc[5L],
    missing_cpg_annotation_records = qc[6L])
  if (qc$invalid_mhap_records > 0) {
    message(
      "Discarded ", qc$invalid_mhap_records, " invalid mHap record(s): ",
      qc$cpg_length_mismatch_records, " CpG/string length mismatch, ",
      qc$invalid_haplotype_string_records, " invalid string, ",
      qc$missing_cpg_annotation_records, " without mapped CpGs."
    )
  }
  seen <- which(!vapply(result$delta_T, is.null, logical(1)))
  deltas <- stats::setNames(lapply(seen, function(i) {
    list(T = result$delta_T[[i]], S = result$delta_S[[i]])
  }), names(positions)[seen])
  list(qc = qc, deltas = deltas, read_files = result$read_files)
}

# Dense per-chromosome loci from a long (chr, pos, sample, S, T) table; a
# repeated (position, sample) keeps its last row.
.dense_loci_from_long <- function(dt, sample_names) {
  result <- list()
  for (chromosome in sort(unique(dt$chr))) {
    rows <- which(dt$chr == chromosome)
    pos <- as.numeric(dt$pos[rows])
    unique_pos <- sort(unique(pos))
    index <- cbind(match(pos, unique_pos), match(dt$sample[rows], sample_names))
    if (anyNA(index)) stop("Locus-table samples must be present in `groups`.", call. = FALSE)
    S <- matrix(0, length(unique_pos), length(sample_names),
                dimnames = list(NULL, sample_names))
    T <- S
    S[index] <- as.numeric(dt$S[rows])
    T[index] <- as.numeric(dt$T[rows])
    result[[chromosome]] <- list(pos = unique_pos, S = S, T = T)
  }
  result
}

# Long locus table (sample, chr, pos, T, S, M) of a store, for diagnostics
# and tests. `candidates = TRUE` returns the screened candidate loci.
.store_locus_table <- function(store, candidates = FALSE) {
  loci <- if (candidates) store$candidates else store$loci
  pieces <- lapply(names(loci), function(chromosome) {
    x <- loci[[chromosome]]
    covered <- which(x$T > 0, arr.ind = TRUE)
    if (!length(covered)) return(NULL)
    if (is.null(dim(covered))) covered <- matrix(covered, ncol = 2L)
    data.table::data.table(sample = colnames(x$T)[covered[, 2L]], chr = chromosome,
      pos = x$pos[covered[, 1L]], T = x$T[covered], S = x$S[covered])
  })
  dt <- data.table::rbindlist(pieces)
  if (!nrow(dt)) {
    return(data.table::data.table(sample = character(), chr = character(),
      pos = numeric(), T = numeric(), S = numeric(), M = numeric()))
  }
  data.table::setorder(dt, chr, pos, sample)
  dt[, M := S / T]
  dt
}

# Region N_sq from the compact read cache: one C sweep per sample/chromosome.
.region_nsq_compact <- function(store, rt, sample_names) {
  N <- matrix(0, nrow(rt), length(sample_names))
  empty <- 0
  for (s in seq_along(sample_names)) {
    paths <- store$read_paths[[sample_names[s]]]
    if (is.null(paths)) next
    key <- store$annotation_by_sample[[sample_names[s]]]
    for (chromosome in intersect(unique(rt$chr), names(paths))) {
      pp <- store$annotations[[key]][[chromosome]]
      rr <- rt[chr == chromosome]
      reads <- .read_compact_reads(paths[[chromosome]])
      result <- .Call(C_region_nsq, reads$lo, reads$hi, reads$a, reads$b,
        reads$count, reads$nd, as.numeric(rr$start), as.numeric(rr$end),
        findInterval(rr$start - 1, pp) + 1L, findInterval(rr$end, pp))
      N[rr$id, s] <- result$N
      empty <- empty + result$empty
    }
  }
  list(N = N, empty = empty)
}

.screen_mscore_store <- function(store, groups, group_levels, min_in_span) {
  # Screen once under observed labels; the full loci stay for region sums.
  store$candidates <- list()
  coverage <- list()
  density <- numeric()
  counts <- c(candidate_loci_before_screen = 0, candidate_loci_retained = 0,
              candidate_loci_dropped_group_coverage = 0)
  for (chromosome in names(store$loci)) {
    x <- store$loci[[chromosome]]
    samples <- colnames(x$T)
    in_group <- function(level) samples %in% names(groups)[groups == level]
    keep <- rowSums(x$T[, in_group(group_levels[1L]), drop = FALSE]) > 0 &
      rowSums(x$T[, in_group(group_levels[2L]), drop = FALSE]) > 0
    dropped <- sum(!keep)
    message("Candidate coverage screen: retained ", sum(keep), " of ",
            length(keep), " CpGs; excluded ", dropped,
            " without positive pooled coverage in both observed groups.")
    counts <- counts + c(length(keep), sum(keep), dropped)
    if (dropped) {
      x <- list(pos = x$pos[keep], S = x$S[keep, , drop = FALSE],
                T = x$T[keep, , drop = FALSE])
    }
    store$candidates[[chromosome]] <- x
    if (length(x$pos)) {
      coverage[[length(coverage) + 1L]] <- rowSums(x$T) / length(groups)
      density <- c(density, min_in_span * (max(x$pos) - min(x$pos) + 1) / length(x$pos))
    }
  }
  all_cov <- unlist(coverage, use.names = FALSE)
  store$cov_q75 <- if (length(all_cov)) as.numeric(stats::quantile(
    all_cov[is.finite(all_cov) & all_cov > 0], .75, type = 7, names = FALSE)) else NA_real_
  store$density_span <- if (length(density)) mean(density) else NA_real_
  store$qc <- c(store$qc, as.list(counts))
  store
}

.smooth_mscore_store <- function(store, groups, group_levels, threshold,
                                 max_gap_smooth, max_gap_region, min_cpgs,
                                 bp_span, min_in_span, check_group_coverage = TRUE,
                                 smoother = "direct") {
  .validate_smoothing_parameters(
    threshold, max_gap_smooth, max_gap_region, min_cpgs, bp_span, min_in_span
  )
  pieces <- list()
  for (chromosome in names(store$candidates)) {
    x <- store$candidates[[chromosome]]
    if (!length(x$pos)) next
    sample_group <- match(unname(groups[colnames(x$T)]), group_levels)
    if (anyNA(sample_group)) {
      stop("Locus-table samples must be present in `groups`.", call. = FALSE)
    }
    locus_stats <- .smoothing_stats_table(chromosome, x, sample_group, store$cov_q75)
    pieces[[length(pieces) + 1L]] <- .smooth_locus_stats(
      locus_stats, threshold, max_gap_smooth, max_gap_region, min_cpgs,
      bp_span, min_in_span, check_group_coverage, store$density_span, smoother)
  }
  pieces <- pieces[lengths(pieces) > 0L]
  if (!length(pieces)) return(GenomicRanges::GRanges())
  frame <- data.table::rbindlist(lapply(pieces, as.data.frame))
  result <- GenomicRanges::GRanges(frame$seqnames, IRanges::IRanges(frame$start, frame$end))
  cols <- setdiff(names(frame), c("seqnames", "start", "end", "width", "strand"))
  S4Vectors::mcols(result) <- S4Vectors::DataFrame(as.data.frame(frame[, cols, with = FALSE]))
  S4Vectors::metadata(result) <- S4Vectors::metadata(pieces[[1L]])
  result
}

# Query finite physical windows, not the full chromosome bounding interval.
# A read-region pair belongs to the window containing max(read_start, region_start).
# This preserves true duplicate input records but prevents repeated window hits.
# Windows bound genomic span, NOT record count; scanTabix(param=...) materializes
# a query. chunk_nrows only bounds parsing/mapping of those returned lines.
.original_region_contributions <- function(source, rt, positions, chunk_nrows,
                                            window_bp = 100000L) {
  if (!isTRUE(all.equal(source$fingerprint,
      file.info(c(source$path, source$index))[, c("size", "mtime"), drop = FALSE]))) {
    stop("Original mHap or index changed during analysis: ", source$path)
  }
  N <- numeric(nrow(rt))
  empty_hits <- list()
  tab <- Rsamtools::TabixFile(source$path, index = source$index)
  open(tab)
  on.exit(close(tab), add = TRUE)
  for (chromosome in intersect(unique(rt$chr), source$chromosomes)) {
    rr <- data.table::copy(rt[chr == chromosome])
    pp <- positions[[chromosome]]
    if (!length(pp)) next
    rr[, `:=`(a_region = findInterval(start - 1, pp) + 1L,
               b_region = findInterval(end, pp))]
    data.table::setkey(rr, chr, start, end)
    window_ids <- sort(unique(unlist(lapply(seq_len(nrow(rr)), function(i)
      seq.int(floor((rr$start[i] - 1) / window_bp),
              floor((rr$end[i] - 1) / window_bp))), use.names = FALSE)))
    for (window_id in window_ids) {
      left <- window_id * window_bp + 1
      right <- min((window_id + 1) * window_bp, .Machine$integer.max - 1)
      # One bp padding accommodates either generic 1-based or BED-style index
      # conventions. Exact physical overlap below always uses mHap coordinates.
      query <- GenomicRanges::GRanges(chromosome,
        IRanges::IRanges(max(1, left - 1), right + 1))
      lines <- Rsamtools::scanTabix(tab, param = query)[[1L]]
      if (!length(lines)) next
      targets <- rr[start <= right & end >= left]
      for (first in seq.int(1L, length(lines), by = chunk_nrows)) {
        raw <- .parse_mhap_lines(lines[seq.int(first, min(length(lines), first + chunk_nrows - 1))])
        raw[, `:=`(sample = source$sample, annotation_key = source$key)]
        # Reuse exact validation/mapping, but skip global methylation state:
        # region S/T already come from the complete locus table. Query QC must
        # not be added again to the ingestion totals.
        reads <- .map_reads_to_cpg(raw, stats::setNames(list(positions), source$key),
                                   compute_global = FALSE)$reads
        if (!nrow(reads)) next
        hits <- data.table::foverlaps(reads, targets,
          by.x = c("chr", "read_start", "read_end"),
          by.y = c("chr", "start", "end"), type = "any", nomatch = NULL)
        if (!nrow(hits)) next
        hits <- hits[pmax(read_start, start) >= left & pmax(read_start, start) <= right]
        hits[, n := pmax(0, pmin(cpg_index_end, b_region) -
                              pmax(cpg_index_start, a_region) + 1)]
        bad <- hits[n == 0, .(id, lo = read_start, hi = read_end,
                              hap = original_hap_string, strand)]
        if (nrow(bad)) {
          bad[, sample := source$sample]
          empty_hits[[length(empty_hits) + 1L]] <- bad
        }
        sums <- hits[n > 0, .(N = sum(count * n^2)), by = id]
        N[sums$id] <- N[sums$id] + sums$N
      }
    }
  }
  list(N = N, empty_hits = empty_hits)
}

.region_matrices_store <- function(store, regions, sample_names) {
  nr <- length(regions)
  dn <- list(if (nr) paste0("region", seq_len(nr)) else character(), sample_names)
  S <- T <- N <- matrix(0, nr, length(sample_names), dimnames = dn)
  rt <- data.table::data.table(chr = as.character(GenomicRanges::seqnames(regions)),
    start = as.numeric(IRanges::start(regions)), end = as.numeric(IRanges::end(regions)),
    id = seq_len(nr))
  # Counts are summed over ALL valid covered loci, not the screened universe.
  for (chromosome in unique(rt$chr)) {
    x <- store$loci[[chromosome]]
    if (is.null(x)) next
    rr <- rt[chr == chromosome]
    sums <- .Call(C_region_sums, x$pos, x$S, x$T, rr$start, rr$end)
    columns <- match(sample_names, colnames(x$T))
    present <- which(!is.na(columns))
    S[rr$id, present] <- sums$S[, columns[present], drop = FALSE]
    T[rr$id, present] <- sums$T[, columns[present], drop = FALSE]
  }
  if (identical(store$read_backend, "compact")) {
    compact <- .region_nsq_compact(store, rt, sample_names)
    N[] <- compact$N
    return(.finish_region_matrices(regions, S, T, N, sample_names,
                                   as.integer(compact$empty)))
  }
  empty_hits <- list()
  for (source in store$read_sources) {
    s <- match(source$sample, sample_names)
    if (is.na(s)) next
    contribution <- .original_region_contributions(source, rt,
      store$annotations[[source$key]], store$chunk_nrows)
    N[, s] <- contribution$N
    empty_hits <- c(empty_hits, contribution$empty_hits)
  }
  ann_key <- NULL
  positions <- NULL
  for (chunk in store$read_chunks) {
    relevant <- rt[chr %in% chunk$spans$chr]
    if (!nrow(relevant)) next
    span_index <- match(relevant$chr, chunk$spans$chr)
    relevant <- relevant[which(relevant$start <= chunk$spans$hi[span_index] &
                                 relevant$end >= chunk$spans$lo[span_index])]
    if (!nrow(relevant)) next
    if (!identical(ann_key, chunk$key)) {
      positions <- store$annotations[[chunk$key]]
      ann_key <- chunk$key
    }
    # One bounding query per chromosome per <=chunk_nrows indexed shard.
    # No record can be returned twice by overlapping queries.
    spans <- relevant[, .(start = min(start), end = max(end)), by = chr]
    query <- GenomicRanges::GRanges(spans$chr, IRanges::IRanges(spans$start, spans$end))
    lines <- unlist(Rsamtools::scanTabix(chunk$path, param = query), use.names = FALSE)
    if (!length(lines)) next
    reads <- data.table::fread(text = paste0(paste(lines, collapse = "\n"), "\n"),
      sep = "\t", header = FALSE,
      col.names = c("chr", "lo", "hi", "a", "b", "count", "global", "hap", "strand"),
      colClasses = c("character", rep("numeric", 5), "logical", "character", "character"),
      showProgress = FALSE)
    for (chromosome in unique(relevant$chr)) {
      pp <- positions[[chromosome]]
      rr <- relevant[chr == chromosome]
      rr[, `:=`(a_region = findInterval(start - 1, pp) + 1L,
                 b_region = findInterval(end, pp))]
      data.table::setkey(rr, chr, start, end)
      hits <- data.table::foverlaps(reads[chr == chromosome], rr,
        by.x = c("chr", "lo", "hi"), by.y = c("chr", "start", "end"),
        type = "any", nomatch = NULL)
      if (!nrow(hits)) next
      hits[, n := pmax(0, pmin(b, b_region) - pmax(a, a_region) + 1)]
      bad <- hits[n == 0, .(id, lo, hi, hap, strand)]
      if (nrow(bad)) {
        bad[, sample := chunk$sample]
        empty_hits[[length(empty_hits) + 1L]] <- bad
      }
      sums <- hits[n > 0, .(value = sum(as.numeric(count) * n^2)), by = id]
      s <- match(chunk$sample, sample_names)
      N[sums$id, s] <- N[sums$id, s] + sums$value
    }
  }
  .finish_region_matrices(regions, S, T, N, sample_names,
    if (length(empty_hits)) nrow(unique(data.table::rbindlist(empty_hits))) else 0L)
}

.finish_region_matrices <- function(regions, S, T, N, sample_names,
                                    reads_without_local_cpg) {
  complete <- rowSums(is.finite(S) & is.finite(T) & is.finite(N) & T > 0 & N > 0) == length(sample_names)
  dropped <- sum(!complete)
  if (dropped) message("  Dropping ", dropped, " region(s) without CpG coverage in every sample.")
  list(regions = regions[complete], S_matrix = S[complete, , drop = FALSE],
    T_matrix = T[complete, , drop = FALSE], N_sq_matrix = N[complete, , drop = FALSE],
    qc = list(reads_without_local_cpg = reads_without_local_cpg,
      regions_dropped_incomplete_coverage = dropped))
}
