#include <R.h>
#include <Rinternals.h>
/* No fused multiply-add contraction (clang contracts a*b + c by default on
 * arm64), so that results are identical across compilers and CPUs. */
#ifdef __clang__
#pragma clang fp contract(off)
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

/*
 * Single-pass mHap reader. For one (optionally gzip/bgzip compressed) mHap
 * file it reproduces the R streaming path exactly:
 *   - record validation and QC counts (as in .map_reads_to_cpg),
 *   - CpG mapping by binary search (as findInterval),
 *   - per-chromosome difference arrays of T and S (as .add_cpg_events),
 *   - optionally the compact read cache (as .consolidate_read_pieces): rows
 *     (lo, hi, a, b, count, nd) per chromosome, sorted by (lo, hi), where nd
 *     counts the distinct (haplotype, strand) records sharing the span.
 * Writing the read cache requires input sorted by chromosome blocks and start;
 * otherwise the scan stops and reports `unsorted` so that R can fall back.
 * Errors are returned as a string after all resources are released.
 */

typedef struct {
    double hi, count;
    int a, b;
    size_t hap_offset;
    int hap_length;
    char strand;
} run_record;

typedef struct {
    /* records sharing the current (chromosome, lo) */
    run_record *records;
    size_t n, capacity;
    char *haps;
    size_t haps_used, haps_capacity;
    /* compact rows of the current chromosome */
    double *lo, *hi, *count, *nd;
    int *a, *b;
    size_t rows, rows_capacity;
} read_buffers;

static const char *sort_haps;

static int compare_run(const void *x, const void *y)
{
    const run_record *p = (const run_record *) x, *q = (const run_record *) y;
    int length, cmp;
    if (p->hi != q->hi) return p->hi < q->hi ? -1 : 1;
    length = p->hap_length < q->hap_length ? p->hap_length : q->hap_length;
    cmp = memcmp(sort_haps + p->hap_offset, sort_haps + q->hap_offset, (size_t) length);
    if (cmp) return cmp;
    if (p->hap_length != q->hap_length) return p->hap_length < q->hap_length ? -1 : 1;
    if (p->strand != q->strand) return p->strand < q->strand ? -1 : 1;
    return 0;
}

static int same_record(const run_record *p, const run_record *q, const char *haps)
{
    return p->hi == q->hi && p->strand == q->strand &&
        p->hap_length == q->hap_length &&
        memcmp(haps + p->hap_offset, haps + q->hap_offset, (size_t) p->hap_length) == 0;
}

static int grow(void **pointer, size_t *capacity, size_t needed, size_t size)
{
    size_t next;
    void *resized;
    if (needed <= *capacity) return 1;
    next = *capacity ? *capacity : 1024;
    while (next < needed) next *= 2;
    resized = realloc(*pointer, next * size);
    if (!resized) return 0;
    *pointer = resized;
    *capacity = next;
    return 1;
}

static int grow_rows(read_buffers *buffers, size_t needed)
{
    size_t capacity = buffers->rows_capacity, next;
    if (needed <= capacity) return 1;
    next = capacity ? capacity : 4096;
    while (next < needed) next *= 2;
#define RESIZE(field, type) do { \
        type *resized = (type *) realloc(buffers->field, next * sizeof(type)); \
        if (!resized) return 0; \
        buffers->field = resized; } while (0)
    RESIZE(lo, double); RESIZE(hi, double); RESIZE(count, double);
    RESIZE(nd, double); RESIZE(a, int); RESIZE(b, int);
#undef RESIZE
    buffers->rows_capacity = next;
    return 1;
}

/* Collapse the current run into compact rows ordered by hi. */
static int flush_run(read_buffers *buffers, double lo)
{
    size_t i = 0, j;
    if (!buffers->n) return 1;
    sort_haps = buffers->haps;
    qsort(buffers->records, buffers->n, sizeof(run_record), compare_run);
    while (i < buffers->n) {
        double count = 0.0, distinct = 0.0;
        for (j = i; j < buffers->n && buffers->records[j].hi == buffers->records[i].hi; j++) {
            count += buffers->records[j].count;
            if (j == i || !same_record(&buffers->records[j], &buffers->records[j - 1], buffers->haps))
                distinct += 1.0;
        }
        if (!grow_rows(buffers, buffers->rows + 1)) return 0;
        buffers->lo[buffers->rows] = lo;
        buffers->hi[buffers->rows] = buffers->records[i].hi;
        buffers->a[buffers->rows] = buffers->records[i].a;
        buffers->b[buffers->rows] = buffers->records[i].b;
        buffers->count[buffers->rows] = count;
        buffers->nd[buffers->rows] = distinct;
        buffers->rows++;
        i = j;
    }
    buffers->n = 0;
    buffers->haps_used = 0;
    return 1;
}

/* Binary layout read by .read_compact_reads(): n, lo, hi, a, b, count, nd. */
static int write_rows(read_buffers *buffers, const char *path)
{
    FILE *file = fopen(path, "wb");
    double n = (double) buffers->rows;
    size_t rows = buffers->rows;
    int ok;
    if (!file) return 0;
    ok = fwrite(&n, sizeof(double), 1, file) == 1 &&
        fwrite(buffers->lo, sizeof(double), rows, file) == rows &&
        fwrite(buffers->hi, sizeof(double), rows, file) == rows &&
        fwrite(buffers->a, sizeof(int), rows, file) == rows &&
        fwrite(buffers->b, sizeof(int), rows, file) == rows &&
        fwrite(buffers->count, sizeof(double), rows, file) == rows &&
        fwrite(buffers->nd, sizeof(double), rows, file) == rows;
    if (fclose(file) != 0) ok = 0;
    buffers->rows = 0;
    return ok;
}

/* Number of positions <= value (positions sorted ascending), as findInterval. */
static R_xlen_t count_at_most(const double *positions, R_xlen_t n, double value)
{
    R_xlen_t left = 0, right = n, mid;
    while (left < right) {
        mid = left + (right - left) / 2;
        if (positions[mid] <= value) left = mid + 1; else right = mid;
    }
    return left;
}

static char *trim(char *start, char *end)
{
    while (start < end && (*start == ' ' || *start == '\t')) start++;
    while (end > start && (end[-1] == ' ' || end[-1] == '\t')) end--;
    *end = '\0';
    return start;
}

static int parse_number(const char *text, double *value)
{
    char *end;
    if (!*text || strcmp(text, "NA") == 0) return 0;
    *value = strtod(text, &end);
    return *end == '\0' && R_FINITE(*value);
}

static void check_interrupt(void *unused) { R_CheckUserInterrupt(); }

static char *copy_string(const char *text)
{
    size_t length = strlen(text) + 1;
    char *copy = (char *) malloc(length);
    if (copy) memcpy(copy, text, length);
    return copy;
}

SEXP C_scan_mhap(SEXP path_sexp, SEXP chr_names, SEXP positions, SEXP read_prefix)
{
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    int write_reads = read_prefix != R_NilValue;
    R_xlen_t n_chr = XLENGTH(chr_names), c;
    gzFile gz = NULL;
    char *line = NULL, *error_text = NULL;
    size_t line_capacity = 1 << 16;
    double qc[6] = {0, 0, 0, 0, 0, 0}; /* total, valid, invalid, mismatch, bad string, missing */
    int unsorted = 0, last_chr = -1, current_chr = -1, n_written = 0;
    double current_lo = 0.0, lines_read = 0.0;
    int *finished = NULL, *written_chr = NULL;
    char **written_paths = NULL;
    read_buffers buffers;
    SEXP deltas_t, deltas_s, result, names, value;
    int n_protect = 0;

    memset(&buffers, 0, sizeof(buffers));
    PROTECT(deltas_t = allocVector(VECSXP, n_chr)); n_protect++;
    PROTECT(deltas_s = allocVector(VECSXP, n_chr)); n_protect++;
    finished = (int *) calloc((size_t) (n_chr > 0 ? n_chr : 1), sizeof(int));
    written_chr = (int *) calloc((size_t) (n_chr > 0 ? n_chr : 1), sizeof(int));
    written_paths = (char **) calloc((size_t) (n_chr > 0 ? n_chr : 1), sizeof(char *));
    line = (char *) malloc(line_capacity);
    if (!finished || !written_chr || !written_paths || !line) {
        error_text = "Out of memory in the mHap reader.";
        goto done;
    }
    gz = gzopen(path, "rb");
    if (!gz) {
        error_text = "Cannot open mHap file.";
        goto done;
    }
    gzbuffer(gz, 1 << 18);

    for (;;) {
        size_t length = 0;
        char *fields[6], *cursor, *end, *chr, *hap, *strand;
        int n_tabs = 0, n_fields, f, hap_ok, basic_ok, global = 0, hap_length;
        double lo, hi, count;
        R_xlen_t first, last, covered, n_pos;
        const double *pos;
        double *delta_t, *delta_s;

        /* Read one full line, growing the buffer for very long records. */
        if (!gzgets(gz, line, (int) line_capacity)) break;
        length = strlen(line);
        while (length > 0 && line[length - 1] != '\n' && !gzeof(gz)) {
            char *resized;
            line_capacity *= 2;
            resized = (char *) realloc(line, line_capacity);
            if (!resized) { error_text = "Out of memory in the mHap reader."; goto done; }
            line = resized;
            if (!gzgets(gz, line + length, (int) (line_capacity - length))) break;
            length += strlen(line + length);
        }
        while (length > 0 && (line[length - 1] == '\n' || line[length - 1] == '\r'))
            line[--length] = '\0';
        if (++lines_read >= 1e6) {
            lines_read = 0;
            if (!R_ToplevelExec(check_interrupt, NULL)) {
                error_text = "Interrupted.";
                goto done;
            }
        }
        /* Blank lines are skipped, as nzchar(trimws(line)) in R. */
        for (cursor = line; *cursor == ' ' || *cursor == '\t' || *cursor == '\r'; cursor++) ;
        if (*cursor == '\0') continue;
        qc[0] += 1;

        for (cursor = line; *cursor; cursor++) if (*cursor == '\t') n_tabs++;
        /* strsplit() drops one trailing empty field. */
        n_fields = n_tabs + 1 - (length > 0 && line[length - 1] == '\t');
        if (n_fields < 6) {
            error_text = "mHap records must have at least six tab-separated fields.";
            goto done;
        }
        cursor = line;
        for (f = 0; f < 6; f++) {
            end = cursor;
            while (*end && *end != '\t') end++;
            {
                int at_end = *end == '\0';
                fields[f] = trim(cursor, end);
                cursor = at_end ? end : end + 1;
            }
        }
        chr = fields[0]; hap = fields[3]; strand = fields[5];
        hap_length = (int) strlen(hap);
        hap_ok = hap_length > 0 && strcmp(hap, "NA") != 0;
        for (f = 0; hap_ok && f < hap_length; f++) {
            if (hap[f] != '0' && hap[f] != '1') hap_ok = 0;
            else if (hap[f] == '1') global = 1;
        }
        if (!hap_ok) qc[4] += 1;
        basic_ok = hap_ok && *chr && strcmp(chr, "NA") != 0 &&
            parse_number(fields[1], &lo) && parse_number(fields[2], &hi) &&
            lo <= hi && parse_number(fields[4], &count) && count > 0 &&
            (strcmp(strand, "+") == 0 || strcmp(strand, "-") == 0);
        if (!basic_ok) continue;

        if (last_chr < 0 || strcmp(CHAR(STRING_ELT(chr_names, last_chr)), chr) != 0) {
            last_chr = -1;
            for (c = 0; c < n_chr; c++) {
                if (strcmp(CHAR(STRING_ELT(chr_names, c)), chr) == 0) { last_chr = (int) c; break; }
            }
        }
        if (last_chr < 0) { qc[5] += 1; continue; }
        pos = REAL(VECTOR_ELT(positions, last_chr));
        n_pos = XLENGTH(VECTOR_ELT(positions, last_chr));
        first = count_at_most(pos, n_pos, lo - 1) + 1;
        last = count_at_most(pos, n_pos, hi);
        covered = last - first + 1;
        if (covered <= 0) { qc[5] += 1; continue; }
        if (covered != hap_length) { qc[3] += 1; continue; }
        qc[1] += 1;

        if (VECTOR_ELT(deltas_t, last_chr) == R_NilValue) {
            SET_VECTOR_ELT(deltas_t, last_chr, allocVector(REALSXP, n_pos + 1));
            SET_VECTOR_ELT(deltas_s, last_chr, allocVector(REALSXP, n_pos + 1));
            memset(REAL(VECTOR_ELT(deltas_t, last_chr)), 0, (size_t) (n_pos + 1) * sizeof(double));
            memset(REAL(VECTOR_ELT(deltas_s, last_chr)), 0, (size_t) (n_pos + 1) * sizeof(double));
        }
        delta_t = REAL(VECTOR_ELT(deltas_t, last_chr));
        delta_s = REAL(VECTOR_ELT(deltas_s, last_chr));
        delta_t[first - 1] += count;
        delta_t[last] -= count;
        if (global) {
            delta_s[first - 1] += count;
            delta_s[last] -= count;
        }

        if (!write_reads) continue;
        if (last_chr != current_chr) {
            if (finished[last_chr]) { unsorted = 1; goto done; }
            if (current_chr >= 0) {
                char file_path[4096];
                if (!flush_run(&buffers, current_lo)) { error_text = "Out of memory in the mHap reader."; goto done; }
                snprintf(file_path, sizeof(file_path), "%s-%d.bin",
                         CHAR(STRING_ELT(read_prefix, 0)), current_chr + 1);
                written_paths[n_written] = copy_string(file_path);
                if (!written_paths[n_written]) { error_text = "Out of memory in the mHap reader."; goto done; }
                written_chr[n_written++] = current_chr;
                if (!write_rows(&buffers, file_path)) { error_text = "Cannot write the read cache."; goto done; }
                finished[current_chr] = 1;
            }
            current_chr = last_chr;
            current_lo = lo;
        } else if (lo < current_lo) {
            unsorted = 1;
            goto done;
        } else if (lo > current_lo) {
            if (!flush_run(&buffers, current_lo)) { error_text = "Out of memory in the mHap reader."; goto done; }
            current_lo = lo;
        }
        if (!grow((void **) &buffers.records, &buffers.capacity, buffers.n + 1, sizeof(run_record)) ||
            !grow((void **) &buffers.haps, &buffers.haps_capacity,
                  buffers.haps_used + (size_t) hap_length, 1)) {
            error_text = "Out of memory in the mHap reader.";
            goto done;
        }
        memcpy(buffers.haps + buffers.haps_used, hap, (size_t) hap_length);
        buffers.records[buffers.n].hi = hi;
        buffers.records[buffers.n].count = count;
        buffers.records[buffers.n].a = (int) first;
        buffers.records[buffers.n].b = (int) last;
        buffers.records[buffers.n].hap_offset = buffers.haps_used;
        buffers.records[buffers.n].hap_length = hap_length;
        buffers.records[buffers.n].strand = strand[0];
        buffers.haps_used += (size_t) hap_length;
        buffers.n++;
    }
    {
        int status;
        const char *message = gzerror(gz, &status);
        if (status != Z_OK && status != Z_STREAM_END) {
            (void) message;
            error_text = "Cannot read the compressed mHap file (corrupt or truncated).";
            goto done;
        }
    }
    if (write_reads && current_chr >= 0) {
        char file_path[4096];
        if (!flush_run(&buffers, current_lo)) { error_text = "Out of memory in the mHap reader."; goto done; }
        snprintf(file_path, sizeof(file_path), "%s-%d.bin",
                 CHAR(STRING_ELT(read_prefix, 0)), current_chr + 1);
        written_paths[n_written] = copy_string(file_path);
                if (!written_paths[n_written]) { error_text = "Out of memory in the mHap reader."; goto done; }
        written_chr[n_written++] = current_chr;
        if (!write_rows(&buffers, file_path)) { error_text = "Cannot write the read cache."; goto done; }
    }

done:
    if (gz) gzclose(gz);
    free(line);
    free(buffers.records); free(buffers.haps);
    free(buffers.lo); free(buffers.hi); free(buffers.count); free(buffers.nd);
    free(buffers.a); free(buffers.b);
    free(finished);

    PROTECT(result = allocVector(VECSXP, 6)); n_protect++;
    PROTECT(names = allocVector(STRSXP, 6)); n_protect++;
    SET_STRING_ELT(names, 0, mkChar("error"));
    SET_STRING_ELT(names, 1, mkChar("unsorted"));
    SET_STRING_ELT(names, 2, mkChar("qc"));
    SET_STRING_ELT(names, 3, mkChar("delta_T"));
    SET_STRING_ELT(names, 4, mkChar("delta_S"));
    SET_STRING_ELT(names, 5, mkChar("read_files"));
    setAttrib(result, R_NamesSymbol, names);
    if (error_text) SET_VECTOR_ELT(result, 0, mkString(error_text));
    SET_VECTOR_ELT(result, 1, ScalarLogical(unsorted));
    PROTECT(value = allocVector(REALSXP, 6)); n_protect++;
    memcpy(REAL(value), qc, sizeof(qc));
    SET_VECTOR_ELT(result, 2, value);
    SET_VECTOR_ELT(result, 3, deltas_t);
    SET_VECTOR_ELT(result, 4, deltas_s);
    {
        /* Paths named by chromosome; R removes them if the scan failed. */
        SEXP files, file_names;
        int i;
        PROTECT(files = allocVector(STRSXP, n_written)); n_protect++;
        PROTECT(file_names = allocVector(STRSXP, n_written)); n_protect++;
        for (i = 0; i < n_written; i++) {
            SET_STRING_ELT(files, i, mkChar(written_paths[i]));
            SET_STRING_ELT(file_names, i, STRING_ELT(chr_names, written_chr[i]));
            free(written_paths[i]);
        }
        setAttrib(files, R_NamesSymbol, file_names);
        SET_VECTOR_ELT(result, 5, files);
    }
    free(written_paths);
    free(written_chr);
    UNPROTECT(n_protect);
    return result;
}
