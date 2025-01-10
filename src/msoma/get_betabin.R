suppressMessages(library(survcomp))
suppressMessages(library(VGAM))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(bbmle))
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(qvalue))
#options(warn=1)

estBetaParams <- function(mu, var) {
  alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta = alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

MODEL_BETABIN_DP_HQ <- function(DATA){
  BB_RESULTS = tryCatch({mle2(ALT_COUNT~dbetabinom.ab(size=DP,shape1,shape2),
                              data=DATA,
                              method="Nelder-Mead",
                              skip.hessian=TRUE,
                              start=list(shape1=mean(DATA$ALT_COUNT),shape2=mean(DATA$REF_COUNT)),#start=list(shape1=1,shape2=round(mean(DATA$DP))),
                              control=list(maxit=1000))},
                        error = function(e) {(estBetaParams(mean(DATA$ALT_FREQ, na.rm = T),var(DATA$ALT_FREQ, na.rm = T)))}) #ALT_FREQ = ALT_COUNT/DP
  ALPHA = ifelse(is.null(coef(BB_RESULTS)[[1]]),BB_RESULTS[[1]], coef(BB_RESULTS)[[1]])
  BETA = ifelse(is.null(coef(BB_RESULTS)[[2]]),BB_RESULTS[[2]], coef(BB_RESULTS)[[2]])

  return (c(ALPHA, BETA))
}

get_alpha_beta <- function(bef_aft,nuclt_change) {
  ref = str_sub(nuclt_change, 1, 1)
  alt = str_sub(nuclt_change, 3, 3)
  refs = c(DNAStringSet(ref), complement(DNAStringSet(ref)))
  nuclt_change_for = paste(DNAStringSet(ref),DNAStringSet(alt),sep='>')
  nuclt_change_rev = paste(complement(DNAStringSet(ref)),complement(DNAStringSet(alt)),sep='>')
  nuclt_changes = c(nuclt_change_for, nuclt_change_rev)
  label = paste0(c(nuclt_change_for, nuclt_change_rev), collapse ="+")
  train_dat = dat %>%
    #filter(BEF_AFT == bef_aft & REF2ALT %in% nuclt_changes & ALT_COUNT > 0)  #AC>0, excluding invariants
    filter(BEF_AFT == bef_aft & #matched context
             ((REF2ALT %in% nuclt_changes & ALT_COUNT > 0) | #AC>0
                 (REF %in% refs & ALT_COUNT == 0)))         #or AC=0
  ## estimate alpha and beta
  fit_bb = MODEL_BETABIN_DP_HQ(train_dat)
  fit_bb_results = data.frame(bef_aft = bef_aft,
                              nuclt_change = nuclt_change,
                              nuclt_change_label = label,
                              nloci = nrow(train_dat),
                              alpha = (fit_bb)[[1]],
                              beta = (fit_bb)[[2]])
  return(fit_bb_results)
}

P_VAL_testacminus0 <- function(DATA, a, b){
  p_val <- pzoibetabinom.ab(DATA$ALT_COUNT-1, DATA$DP, a, b, lower.tail=F)
  p_val[p_val < 0] <- 0
  return(p_val)
}

P_VAL_testacminus1 <- function(DATA, a, b){
  p_val <- pzoibetabinom.ab(DATA$ALT_COUNT-2, DATA$DP-1, a, b, lower.tail=F)
  p_val[p_val < 0] <- 0
  return(p_val)
}

P_VAL_testacminus2 <- function(DATA, a, b){
  p_val <- pzoibetabinom.ab(DATA$ALT_COUNT-3, DATA$DP-2, a, b, lower.tail=F)
  p_val[p_val < 0] <- 0
  return(p_val)
}

P_VAL_testacminus3 <- function(DATA, a, b){
  p_val <- pzoibetabinom.ab(DATA$ALT_COUNT-4, DATA$DP-3, a, b, lower.tail=F)
  p_val[p_val < 0] <- 0
  return(p_val)
}

get_adj_p <- function(bef_aft,nuclt_change) {
  ref = str_sub(nuclt_change, 1, 1)
  alt = str_sub(nuclt_change, 3, 3)
  refs = c(DNAStringSet(ref), complement(DNAStringSet(ref)))
  nuclt_change_for = paste(DNAStringSet(ref),DNAStringSet(alt),sep='>')
  nuclt_change_rev = paste(complement(DNAStringSet(ref)),complement(DNAStringSet(alt)),sep='>')
  nuclt_changes = c(nuclt_change_for, nuclt_change_rev)
  label = paste0(nuclt_changes,collapse ="+")

  alpha_beta = train_alpha_beta %>%
    filter(bef_aft == bef_aft & nuclt_change == nuclt_change) #minDP == args$mindp & maxDP == args$maxdp & maxAltFreq == args$maxaltfreq &
  alpha = unlist(alpha_beta$alpha[1])
  beta = unlist(alpha_beta$beta[1])

  calls_dat = dat %>%
    filter(BEF_AFT == bef_aft & REF2ALT %in% nuclt_changes)

  #RB dropping rows with null values (maybe don't do this?) had to do it on small test datasets
  calls_dat = na.omit(calls_dat)

  #RB if the table is too small to adjust pvalues on
  #RB otherwise we get errors trying to get p-values
  if(nrow(calls_dat) <= 1){
      return(calls_dat)
  }

  calls_dat$P_VAL <- P_VAL_testacminus0(calls_dat, alpha, beta)
  calls_dat$P_VAL_adj <- p.adjust(calls_dat$P_VAL, method = "bonferroni")
  calls_dat$P_VAL_acm0_BHadj <- p.adjust(calls_dat$P_VAL, method = "BH")
  #qobj_acm0 <- qvalue(p = calls_dat$P_VAL_acm0)
  #calls_dat$qvalues_acm0 <- qobj_acm0$qvalues
  #calls_dat$pi0_acm0 <- qobj_acm0$pi0
  #calls_dat$lfdr_acm0 <- qobj_acm0$lfdr

  calls_dat$P_VAL_acm1 <- P_VAL_testacminus1(calls_dat, alpha, beta)
  calls_dat$P_VAL_acm1_BFadj <- p.adjust(calls_dat$P_VAL_acm1, method = "bonferroni")
  calls_dat$P_VAL_acm1_BHadj <- p.adjust(calls_dat$P_VAL_acm1, method = "BH")
  qobj_acm1 <- qvalue(p = calls_dat$P_VAL_acm1)
  calls_dat$qvalues_acm1 <- qobj_acm1$qvalues
  calls_dat$pi0_acm1 <- qobj_acm1$pi0
  calls_dat$lfdr_acm1 <- qobj_acm1$lfdr

  calls_dat$P_VAL_acm2 <- P_VAL_testacminus2(calls_dat, alpha, beta)
  calls_dat$P_VAL_acm2_BFadj <- p.adjust(calls_dat$P_VAL_acm2, method = "bonferroni")
  calls_dat$P_VAL_acm2_BHadj <- p.adjust(calls_dat$P_VAL_acm2, method = "BH")
  qobj_acm2 <- qvalue(p = calls_dat$P_VAL_acm2)
  calls_dat$qvalues_acm2 <- qobj_acm2$qvalues
  calls_dat$pi0_acm2 <- qobj_acm2$pi0
  calls_dat$lfdr_acm2 <- qobj_acm2$lfdr

  calls_dat$P_VAL_acm3 <- P_VAL_testacminus3(calls_dat, alpha, beta)
  calls_dat$P_VAL_acm3_BFadj <- p.adjust(calls_dat$P_VAL_acm3, method = "bonferroni")
  calls_dat$P_VAL_acm3_BHadj <- p.adjust(calls_dat$P_VAL_acm3, method = "BH")
  qobj_acm3 <- qvalue(p = calls_dat$P_VAL_acm3)
  calls_dat$qvalues_acm3 <- qobj_acm3$qvalues
  calls_dat$pi0_acm3 <- qobj_acm3$pi0
  calls_dat$lfdr_acm3 <- qobj_acm3$lfdr

  return(calls_dat)
}

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-input", "--input", type="character", help="Input to trian error rate distribution", nargs=1, required=TRUE)
parser$add_argument("-output", "--output", type="character", help="Output with Beta Binomial P values", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-ab_file", "--ab_file", type="character", help="Files with trained alpha and beta", nargs=1, required=TRUE)
parser$add_argument("-mindp", "--mindp", type="integer", help="minimum DP to filter the call sets", required=TRUE)
#maxdp is disabled; set to be sample-specific, 5*median_dp
#parser$add_argument("-maxdp", "--maxdp", type="integer", help="maximum DP to filter the call sets", required=TRUE)

# reading parameters
# args$input == The counts file
# args$output == Output path for p-values
# args$ab_file == Output path for a/b trained params
# args$mindp == 10
args = parser$parse_args()

dat = fread(args$input) %>%
  as_tibble() %>%
  mutate(REF_COUNT = str_replace_all(REF_COUNT, "[.]", "0"),
         ALT_COUNT = str_replace_all(ALT_COUNT, "[.]", "0"),
         REF_COUNT = as.integer(REF_COUNT),
         ALT_COUNT = as.integer(ALT_COUNT),
         DP = REF_COUNT + ALT_COUNT)

maxdp = 5*median(dat$DP)
dat = dat %>%
  filter(DP >= args$mindp & DP <= maxdp) %>%
  mutate(ALT_FREQ = ALT_COUNT/DP,
         BEF_AFT = ifelse(REF %in% c("C","T"), paste(BEFORE,AFTER,sep=""), paste(complement(DNAStringSet(AFTER)),complement(DNAStringSet(BEFORE)),sep="")),
         REF2ALT = paste(REF,ALT,sep='>'))

nuclt_changes_6 = c("T>G","T>A","T>C","C>G","C>T","C>A")
bases = c("A","T","C","G")
params = expand.grid(bases, bases, nuclt_changes_6) %>%
  unite("bef_aft", Var1, Var2,sep = "") %>%
  dplyr::rename(nuclt_change = Var3)

train_alpha_beta <- params %>%
  pmap(get_alpha_beta) %>%
  bind_rows() %>%
  as_tibble() %>%
  mutate(minDP = args$mindp,
         maxDP = maxdp)


if (grepl('gz$', args$ab_file)) {
  write.table(train_alpha_beta, gzfile(args$ab_file), row.names=FALSE, quote=FALSE, sep="\t")
} else {
  write.table(train_alpha_beta, args$ab_file, row.names=FALSE, quote=FALSE, sep="\t")
}


calls_wP = params %>%
  pmap(get_adj_p) %>%
  bind_rows()

## Fisher strand bias
#calls_wP$rawFISHER = apply(calls_wP,1,function(x) fisher.test(matrix(c(as.numeric(x["REF_FWD"]), as.numeric(x["ALT_FWD"]), as.numeric(x["REF_REV"]), as.numeric(x["ALT_REV"])), nrow = 2))[[1]])
#calls_wP$adjFISHER = p.adjust(calls_wP$rawFISHER, method = "bonferroni")

if (grepl('gz$', args$output)) {
  write.table(calls_wP, gzfile(args$output), row.names=FALSE, quote=FALSE, sep="\t")
} else {
  write.table(calls_wP, args$output, row.names=FALSE, quote=FALSE, sep="\t")
}
