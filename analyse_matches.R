library(ggplot2)
library(dplyr)
library(tibble)
library(vroom)

matches_df <- vroom::vroom("test_csv_real.csv", col_types = cols(comp = col_character(), .default = col_double())) %>%
  column_to_rownames("comp")

scale_lb <- seq(0.5, 35.5, by = 3)
cur_pow <- 0.1
power_increment <- 0.05
ls <- length(scale_lb)
while (scale_lb[length(scale_lb)] < max(as.integer(colnames(matches_df)))) {
  scale_lb <- c(scale_lb, scale_lb[ls] * 10^cur_pow)
  cur_pow <- cur_pow + power_increment
}

all_matches <- colSums(matches_df, na.rm = TRUE)
all_matches_df <- tibble("match_length" = as.integer(names(all_matches)), "freq" = all_matches) %>%
  arrange(match_length)


binned_match_df <- tibble(x = scale_lb, y = rep(0, length.out = length(scale_lb)))
i <- 1
j <- 1

while(i <= nrow(all_matches_df)) {
  if(all_matches_df[i, "match_length"] <= binned_match_df[j, "x"]) {
    binned_match_df[j, "y"] <- binned_match_df[j, "y"] + all_matches_df[i, "freq"]
    i <- i + 1
  } else {
    j <- j + 1
    if(j > nrow(binned_match_df))
      break
  }
}

# il faut normaliser
for(r_i in seq_len(nrow(binned_match_df) - 1)) {
  len_bin <- binned_match_df[r_i + 1, "x"] - binned_match_df[r_i, "x"]
  binned_match_df[r_i + 1, "y"] <- binned_match_df[r_i + 1, "y"] / len_bin
}

# ok ça a un peu une tête de double loi puissance imo
ggplot(all_matches_df, aes(x = match_length, y = freq)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

ggsave("plot_fig2_2.png")

# optimization -----------------------------------------------------------------

# some defs : (fuck that but no ... in metaops)
mus <- 5e-9
muc <- 6e-11
d <- 0.55
r <- binned_match_df$x[2:nrow(binned_match_df)]
mE <- binned_match_df$y[2:nrow(binned_match_df)]
L0 <- 4903888.5 * 10^4
L0_fit <- FALSE


theoretical_MLD <- function(par1){
  dr <- 0.1
  rm <- r - dr
  rp <- r + dr
  tau <- as.numeric(10^par1[1])
  rho <- as.numeric(10^par1[2])
  if (length(par1) >2 ){
    print(par1)
    print(par1[3])
    if ( is.na(par1[3]) == F ){
      L0 = 10^par1[3]
    }
  }
  mua <- min(d/tau,mus)

  mc <- 2*((1 + r*mua*tau)/exp(r*mua*tau) - (1 + r*muc*tau)/exp(r*muc*tau))/(r^2*(muc^2 - mus^2)*tau^2)
  mcm <- 2*((1 + rm*mua*tau)/exp(rm*mua*tau) - (1 + rm*muc*tau)/exp(rm*muc*tau))/(rm^2*(muc^2 - mus^2)*tau^2)
  mcp <- 2*((1 + rp*mua*tau)/exp(rp*mua*tau) - (1 + rp*muc*tau)/exp(rp*muc*tau))/(rp^2*(muc^2 - mus^2)*tau^2)


  mc <- L0*(mcm + mcp - 2*mc)/dr^2
  mc[is.na(mc)] <- 0

  if (tau<d/mus)
  {
    mh <- (2*(-exp(-(r*muc*tau)) + exp(-(r*mus*tau)) + r*(-muc + mus)*tau))/(r^2*(-muc^2 + mus^2)*tau)
    mhm <- (2*(-exp(-(rm*muc*tau)) + exp(-(rm*mus*tau)) + rm*(-muc + mus)*tau))/(rm^2*(-muc^2 + mus^2)*tau)
    mhp <- (2*(-exp(-(rp*muc*tau)) + exp(-(rp*mus*tau)) + rp*(-muc + mus)*tau))/(rp^2*(-muc^2 + mus^2)*tau)
  }else
  {
    mh <- (-2*(-exp(-(r*muc*tau)) + r*(-muc + mus)*tau + (1 + r*(d - mus*tau))/exp(r*d)))/(r^2*(muc^2 - mus^2)*tau)
    mhm <- (-2*(-exp(-(rm*muc*tau)) + rm*(-muc + mus)*tau + (1 + rm*(d - mus*tau))/exp(rm*d)))/(rm^2*(muc^2 - mus^2)*tau)
    mhp <- (-2*(-exp(-(rp*muc*tau)) + rp*(-muc + mus)*tau + (1 + rp*(d - mus*tau))/exp(rp*d)))/(rp^2*(muc^2 - mus^2)*tau)
  }
  mh <- L0*rho*(mhm + mhp - 2*mh)/dr^2

  return(list("mh" = mh, "mc" = mc))
}



LLlocal <- function(par1)
{
  mCalc <- theoretical_MLD(par1)

  mT <- mCalc$mc + mCalc$mh
  return(mean( ((mT-mE)/(mE + mT))^2))
}

library(metaheuristicOpt)

if (L0_fit == "yes"){
res <- metaOpt(LLlocal, optimType = "MIN", algorithm = c("HS"), 3,
               rangeVar=t(as.matrix(data.frame(lower=log10(c(1e6,1e-13,1e5)),
                                               upper = log10(c(5e9,1.01e-8,1e7))))),
               control = list(numPopulation=10000,maxIter=10), seed = 3)
}else{
  res <- metaOpt(LLlocal, optimType = "MIN", algorithm = c("HS"), 2,
                 rangeVar=t(as.matrix(data.frame(lower=log10(c(1e6,1e-13)),
                                                 upper = log10(c(5e9,1.01e-8))))),
                 control = list(numPopulation=10000,maxIter=10), seed = 3)
}

print(paste("The estimated value of tau is :", 10^res$result[1], "\n", "The estimated value of rho is :", 10^res$result[2]))
# C'est le bon rho mais loin d'être le bon tau on dirait
#  affreux usage de variables globales
mc_fun <- function(match_length) {
  r <<- match_length
  return(theoretical_MLD(res$result)$mc)
}
mh_fun <- function(match_length) {
  r <<- match_length
  return(theoretical_MLD(res$result)$mh)
}
sum_equa <- function(match_length) {
  return(mh_fun(match_length) + mc_fun(match_length))
}
th_fun <- function(match_length) {
  r <<- match_length
  return((theoretical_MLD(c(8, -10))$mh + theoretical_MLD(c(8, -10))$mc) * 10^4)
}

r4 <- function(match_length) {
  r <<- match_length
  return(mE[1]/match_length^4)
}

# Final representation

ggplot(binned_match_df %>% filter(x != 0.5), aes(x = x, y = y)) +
  geom_point() +
  geom_function(fun = mc_fun, color = "blue") +
  geom_function(fun = mh_fun, color = "red") +
  geom_function(fun = sum_equa, color = "black") +
  geom_function(fun = r4, color = "green") +
  geom_function(fun = th_fun, color = "purple") +
  scale_x_log10() +
  scale_y_log10()


res_opt <- optim(c(log10(10^8), log10(10^-10)),
                 LLlocal,
                 method = "L-BFGS-B",
                 lower = log10(c(1e6,1e-13)),
                 upper = log10(c(5e9,1.01e-8)),
                 )
