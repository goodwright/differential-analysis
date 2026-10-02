# Validation is deliberately separate from fitting so design errors fail before analysis.
read_table <- function(path) {
    first <- readLines(path, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
    read.table(path, header = TRUE, sep = sep, quote = '"', comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE)
}
write_tsv <- function(x, path) {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
parse_names <- function(x) {
    if (is.null(x) || !nzchar(x)) return(character())
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}
validated_design <- function(samples, formula_text, numeric_covariates = "") {
    assert(!anyDuplicated(names(samples)), "Samplesheet has duplicate column names")
    assert("sample_id" %in% names(samples), "Samplesheet requires sample_id")
    assert(!anyNA(samples$sample_id) && all(nzchar(samples$sample_id)) &&
           !anyDuplicated(samples$sample_id), "Sample IDs must be present and unique")
    parsed <- tryCatch(str2lang(formula_text), error = function(e) stop("Invalid design formula"))
    allowed <- c("~", "+", "-", "*", ":", "(")
    validate_ast <- function(x) {
        if (is.symbol(x)) {
            assert(as.character(x) %in% names(samples), paste("Unknown model column:", x))
            assert(as.character(x) != "sample_id", "sample_id is an identifier, not a model factor; provide a separate pair column")
        } else if (is.numeric(x)) {
            assert(length(x) == 1 && x %in% c(0, 1), "Only 0 and 1 are allowed as formula constants")
        } else if (is.call(x)) {
            assert(as.character(x[[1]]) %in% allowed, "Design formula permits columns and ~ + - * : () only; functions are not allowed")
            lapply(as.list(x)[-1], validate_ast)
        } else stop("Unsupported design expression", call. = FALSE)
        invisible(NULL)
    }
    assert(is.call(parsed) && identical(parsed[[1]], as.name("~")) && length(parsed) == 2,
           "Design must be a one-sided formula, e.g. ~ pair + poly + poly:ko1")
    validate_ast(parsed)
    form <- as.formula(parsed, env = baseenv())
    variables <- all.vars(form)
    numeric_names <- parse_names(numeric_covariates)
    assert(all(numeric_names %in% variables), "Numeric covariates must be columns used in the design")
    for (name in variables) {
        x <- samples[[name]]
        assert(!anyNA(x) && all(nzchar(as.character(x))), paste("Missing model values:", name))
        if (name %in% numeric_names) {
            value <- suppressWarnings(as.numeric(x))
            assert(all(is.finite(value)), paste("Non-numeric/nonfinite covariate:", name))
            samples[[name]] <- value
        } else samples[[name]] <- factor(x, levels = sort(unique(as.character(x))))
    }
    design <- model.matrix(form, samples, na.action = na.fail)
    rownames(design) <- samples$sample_id
    assert(all(is.finite(design)), "Nonfinite design matrix")
    rank <- qr(design)$rank
    assert(rank == ncol(design), "Rank-deficient design: remove confounded terms (pair already encodes genotype). No coefficients were dropped automatically.")
    assert(nrow(design) - rank >= 2, "Design requires at least two residual degrees of freedom")
    list(matrix = design, samples = samples, formula = formula_text,
         residual_df = nrow(design) - rank,
         factor_levels = lapply(samples[setdiff(variables, numeric_names)], levels))
}
validated_contrasts <- function(table, coefficients) {
    assert(all(c("contrast", "coefficient", "weight") %in% names(table)),
           "Contrast table requires contrast, coefficient, weight columns")
    assert(nrow(table) > 0 && !anyNA(table[c("contrast", "coefficient", "weight")]), "Empty/missing contrast values")
    assert(all(grepl("^[A-Za-z][A-Za-z0-9_]*$", table$contrast)), "Contrast names must be safe identifiers: letters, digits, underscores")
    assert(!anyDuplicated(table[c("contrast", "coefficient")]), "Duplicate contrast/coefficient rows")
    assert(all(table$coefficient %in% coefficients), paste("Unknown contrast coefficient; available:", paste(coefficients, collapse = ", ")))
    weights <- suppressWarnings(as.numeric(table$weight))
    assert(all(is.finite(weights)), "Contrast weights must be finite numbers")
    names <- unique(table$contrast)
    mat <- matrix(0, length(coefficients), length(names), dimnames = list(coefficients, names))
    mat[cbind(match(table$coefficient, coefficients), match(table$contrast, names))] <- weights
    assert(all(colSums(abs(mat)) > 0), "A contrast cannot be all zero")
    mat
}
read_gmt <- function(path) {
    lines <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
    assert(length(lines) > 0 && all(lengths(lines) >= 3), "GMT requires name, description and gene IDs")
    names <- vapply(lines, `[`, "", 1)
    assert(all(nzchar(names)) && !anyDuplicated(names), "GMT names must be nonempty and unique")
    sets <- lapply(lines, function(x) unique(x[-c(1, 2)]))
    assert(all(vapply(sets, function(x) all(nzchar(x)), TRUE)), "Empty GMT gene ID")
    setNames(sets, names)
}
