#!/usr/bin/env Rscript
# Synthetic data only. Tests model identifiability, contrast direction and nuisance separation.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "lib/advanced_model/helpers.R"))
library(jsonlite)
out <- if (length(args) >= 2) args[2] else tempfile("advanced-model-tests-")
dir.create(out, recursive = TRUE, showWarnings = FALSE); out <- normalizePath(out)
set.seed(48173)
s <- expand.grid(fraction = c("input", "polysome"), rep = 1:2, line = c("WT", "KO1", "KO9"), stringsAsFactors = FALSE)
s$sample_id <- paste(s$line, s$rep, s$fraction, sep = "_"); s$condition <- s$line
s$pair <- paste(s$line, s$rep, sep = "_");s$poly <- as.integer(s$fraction == "polysome")
s$ko1 <- as.integer(s$line == "KO1");s$ko9 <- as.integer(s$line == "KO9")
s <- s[,c("sample_id","condition","pair","fraction","poly","ko1","ko9")]
formula <- "~ 0 + pair + poly + poly:ko1 + poly:ko9"
model <- validated_design(s, formula, "poly,ko1,ko9")
stopifnot(model$residual_df == 3)
ct <- data.frame(contrast = c("average","average","clone1","clone9","reverse","reverse"),
 coefficient = c("poly:ko1","poly:ko9","poly:ko1","poly:ko9","poly:ko1","poly:ko9"),
 weight = c(.5,.5,1,1,-.5,-.5))
c <- validated_contrasts(ct, colnames(model$matrix))
stopifnot(all(c[,"reverse"] == -c[,"average"]))
checks <- 0
must_fail <- function(expr, pattern) {
    e <- tryCatch({force(expr);NULL}, error = identity)
    stopifnot(inherits(e,"error"), grepl(pattern, conditionMessage(e), ignore.case = TRUE))
    checks <<- checks + 1
}
must_fail(validated_design(s,"~ 0 + pair + condition + poly","poly"),"Rank-deficient")
must_fail(validated_design(s,"~ system('touch should-not-exist')",""),"functions are not allowed")
must_fail(validated_design(s,"~ pair + missing",""),"Unknown model column")
s_bad <- s;s_bad$poly[1] <- NA
must_fail(validated_design(s_bad,formula,"poly,ko1,ko9"),"Missing model values")
s_bad <- s;s_bad$sample_id[1] <- s_bad$sample_id[2]
must_fail(validated_design(s_bad,formula,"poly,ko1,ko9"),"unique")
ct_bad <- ct;ct_bad$coefficient[1] <- "not_a_coefficient"
must_fail(validated_contrasts(ct_bad,colnames(model$matrix)),"Unknown contrast coefficient")
ct_bad <- ct;ct_bad$weight[] <- 0
must_fail(validated_contrasts(ct_bad,colnames(model$matrix)),"all zero")
must_fail(validated_contrasts(rbind(ct,ct[1,]),colnames(model$matrix)),"Duplicate")
s_bad <- s;s_bad$poly <- as.character(s_bad$poly);s_bad$poly[1] <- "bad"
must_fail(validated_design(s_bad,formula,"poly,ko1,ko9"),"Non-numeric")
s_bad <- s;s_bad$unique_pair <- seq_len(nrow(s_bad))
must_fail(validated_design(s_bad,"~ 0 + unique_pair",""),"residual degrees")
# Paired data with an input-RNA nuisance effect in a different module.
ng <- 600; mu <- matrix(exp(runif(ng,log(300),log(1800))),ng,nrow(s))
mu[1:30,] <- sweep(mu[1:30,],2,2^(-1.2*s$poly*s$ko1 - 1.0*s$poly*s$ko9),"*")
mu[31:60,] <- sweep(mu[31:60,],2,2^(-1.5*(s$ko1+s$ko9)),"*")
mu <- sweep(mu,2,rep(c(.7,1.3),6),"*")
counts <- matrix(rnbinom(length(mu),mu=as.vector(mu),size=80),ng,nrow(s));colnames(counts)<-s$sample_id
raw <- data.frame(gene_id=paste0("g",seq_len(ng)),counts,check.names=FALSE)
write_tsv(raw,file.path(out,"counts.tsv"));write.csv(s,file.path(out,"samples.csv"),row.names=FALSE)
write.csv(ct,file.path(out,"contrasts.csv"),row.names=FALSE)
writeLines(c(paste(c("translation_module","synthetic",paste0("g",1:30)),collapse="\t"),paste(c("RNA_only_module","synthetic",paste0("g",31:60)),collapse="\t")),file.path(out,"modules.gmt"))
writeLines(raw$gene_id,file.path(out,"universe.txt"))
opt <- list(design_formula=formula,numeric_covariates="poly,ko1,ko9",seed=20261002,
 count_matrix_type="raw_counts",filter_column="fraction",filter_level="input",min_cpm=1,min_samples=2,
 gene_set_id_column="gene_id",use_gene_sets=TRUE,use_universe=TRUE,min_set_size=10,max_set_size=500,set_adjust_method="BH")
write_json(opt,file.path(out,"options.json"),auto_unbox=TRUE,pretty=TRUE)
run <- function(destination, counts_file="counts.tsv", options="options.json") {
 dir.create(file.path(out,destination),showWarnings=FALSE)
 old <- setwd(file.path(out,destination));on.exit(setwd(old))
 paths <- c(file.path(root,"lib/advanced_model/run.R"),file.path(out,c(options,counts_file,"samples.csv","contrasts.csv","modules.gmt","universe.txt")))
 result <- system2(file.path(R.home("bin"),"Rscript"),shQuote(paths),stdout="run.log",stderr="run.log")
 if (result != 0) stop(paste(readLines("run.log"),collapse="\n"))
 read_table("camera.results.tsv")
}
result <- run("result")
get <- function(set,contrast) result[result$gene_set==set & result$contrast==contrast,,drop=FALSE]
a <- get("translation_module","average");b <- get("translation_module","reverse");n <- get("RNA_only_module","average")
stopifnot(a$excess_log2_change < -.7,a$PValue < .01,a$Direction == "Down")
stopifnot(abs(n$excess_log2_change) < .3)
stopifnot(abs(a$excess_log2_change+b$excess_log2_change)<1e-10,b$Direction=="Up",b$PValue<.01)
gene_a <- read_table(file.path(out,"result/average.limma.results.tsv"))
gene_b <- read_table(file.path(out,"result/reverse.limma.results.tsv"))
stopifnot(max(abs(gene_a$P.Value-gene_b$P.Value))<1e-10)
stopifnot(isTRUE(all.equal(result$adjusted_p,p.adjust(result$PValue,"BH"))))
# Count-column order must not change the result.
write_tsv(raw[,c(1,rev(seq.int(2,ncol(raw))))],file.path(out,"shuffled.tsv"))
shuffled <- run("shuffled_result","shuffled.tsv")
stopifnot(isTRUE(all.equal(result,shuffled,tolerance=1e-10)))
# Gene-set coverage failure is explicit rather than a spurious favorable test.
coverage_opts <- opt;coverage_opts$min_set_size <- 40
write_json(coverage_opts,file.path(out,"coverage-options.json"),auto_unbox=TRUE)
dir.create(file.path(out,"coverage_result"),showWarnings=FALSE)
old <- setwd(file.path(out,"coverage_result"))
paths <- c(file.path(root,"lib/advanced_model/run.R"),file.path(out,c("coverage-options.json","counts.tsv","samples.csv","contrasts.csv","modules.gmt","universe.txt")))
status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(paths),stdout="run.log",stderr="run.log")
stopifnot(status==0,file.exists("camera.not_tested.tsv"),!file.exists("camera.results.tsv"));setwd(old)
write_json(list(status="passed",validation_rejection_cases=checks,synthetic_average_excess=a$excess_log2_change,
 synthetic_average_p=a$PValue,RNA_only_excess=n$excess_log2_change,
 verified=c("paired interaction recovers simulated change","total-RNA-only nuisance separated","reverse contrast flips direction and preserves gene-level p-values","count columns reordered safely","module correction across all set/contrast tests","insufficient module coverage recorded")),file.path(out,"test-summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat(readLines(file.path(out,"test-summary.json")),sep="\n")
