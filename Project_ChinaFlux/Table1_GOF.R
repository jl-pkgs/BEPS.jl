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

write_list2xlsx(tbl, "Table1_GOF_V3.xlsx", show = TRUE)
