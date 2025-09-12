## KS test wrapper
options(box.path=file.path(Sys.getenv("CBMGIT"), "MLscripts"))
box::use(KS = R/ks_genescore)

#dat <- read.delim("~/Desktop/Book7.txt") |>
dat <- readxl::read_xlsx("~/Desktop/Worksheet in 2024_10_21_smonti_ziwhuang_hypeR-GEM.v2.xlsx") |>
  dplyr::select(Pvalue) |>
  dplyr::filter(Pvalue != "N/A") |>
  tidyr::separate_wider_delim(Pvalue, delim = " (", names = c("bi", "uni")) |>
  dplyr::mutate(uni = stringr::str_remove(uni, "\\)")) |>
  dplyr::mutate(across(everything(), as.numeric))

ks.test(x = dat$bi, y = "punif", alternative = "greater")
ks.test(x = unique(dat$uni), y = "punif", alternative = "greater")
.21
