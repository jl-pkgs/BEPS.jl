# %% GPP & ET
gof_prev <- map(lst, \(l) l$gof$Flux) %>% melt_list("site")
gof_optim <- map(lst, \(l) l$gof_opt$Flux) %>% melt_list("site")

gof <- list("Original" = gof_prev, "Optimized" = gof_optim) %>% melt_list("type")
dat <- melt(gof, c("type", "site", "var", "n_valid"), variable.name = "index") %>%
  dcast(site + var + index ~ type, value.var = "value") %>%
  merge(st, by = "site")

Table1 <- function(VAR = "GPP", INDEX = "NSE") {
  mutate(dat, diff = Optimized - Original) %>%
    relocate(I, site, var, index, Original) %>%
    arrange(I) %>%
    subset(var == VAR & index == INDEX) %>%
    select(-label, -index) %>%
    dt_round(3)
}

tbl <- list(
  NSE = cbind(Table1("GPP", "NSE"), Table1("ET", "NSE")[, -(1:2)]),
  KGE = cbind(Table1("GPP", "KGE"), Table1("ET", "KGE")[, -(1:2)]),
  R2 = cbind(Table1("GPP", "R2"), Table1("ET", "R2")[, -(1:2)])
)
write_list2xlsx(tbl, "Table1_GOF_V4.xlsx", show = TRUE)

# %% SM
# l <- lst[[1]]
gof_prev <- map(lst, \(l) check_var(l$gof$SM)) %>% melt_list("site")
gof_optim <- map(lst, \(l) check_var(l$gof_opt$SM)) %>% melt_list("site")
gof <- list("Original" = gof_prev, "Optimized" = gof_optim) %>% melt_list("type")

dat <- melt(gof, c("type", "site", "var", "n_valid"), variable.name = "index") %>%
  dcast(site + var + index ~ type, value.var = "value") %>%
  merge(st, by = "site")

Table1 <- function(INDEX = "NSE") {
  mutate(dat, diff = Optimized - Original) %>%
    relocate(I, site, var, index, Original) %>%
    arrange(I) %>%
    subset(index == INDEX) %>%
    select(-label, -index) %>%
    dt_round(3)
}

tbl <- list(
  NSE = Table1("NSE"),
  KGE = Table1("KGE"),
  R2 = Table1("R2")
)
# write_list2xlsx(tbl, "Table1_SM_GOF.xlsx", show = TRUE)

# %% TS
# l <- lst[[1]]
gof_prev <- map(lst, \(l) check_var(l$gof$TS)) %>% melt_list("site")
gof_optim <- map(lst, \(l) check_var(l$gof_opt$TS)) %>% melt_list("site")
gof <- list("Original" = gof_prev, "Optimized" = gof_optim) %>% melt_list("type")

dat <- melt(gof, c("type", "site", "var", "n_valid"), variable.name = "index") %>%
  dcast(site + var + index ~ type, value.var = "value") %>%
  merge(st, by = "site")

Table1 <- function(INDEX) {
  mutate(dat, diff = Optimized - Original) %>%
    relocate(I, site, var, index, Original) %>%
    arrange(I) %>%
    subset(index == INDEX) %>%
    select(-label, -index) %>%
    dt_round(3)
}

tbl <- list(
  NSE = Table1("NSE"),
  KGE = Table1("KGE"),
  R2 = Table1("R2")
)
write_list2xlsx(tbl, "Table1_TS_GOF.xlsx", show = TRUE)
