library(data.table)

read_hmd_wide <- function(path, value_name) {
  dt <- fread(path, skip = 2, na.strings = ".")
  melt(dt, id.vars = c("Year", "Age"),
       measure.vars = c("Female", "Male", "Total"),
       variable.name = "sex", value.name = value_name)
}

# 1. Reference standard raw input: Sweden 2000-2019 (1x1)
swe_mx <- read_hmd_wide("data-raw/hmd/SWE/Mx_1x1.txt", "mx")[Year %in% 2000:2019]
swe_mx[, age := as.integer(sub("\\+", "", Age))]
swe_mx[, x := fifelse(age == 0L, 0.2, fifelse(age == 110L, 111.5, age + 0.5))]
std_input <- swe_mx[is.finite(mx) & mx > 0, .(year = Year, sex = as.character(sex), x, mx)]

# 2. Target country-years
targets <- list(
  list(code = "USA", year = 2015L),
  list(code = "RUS", year = 2000L),
  list(code = "JPN", year = 1970L),
  list(code = "ISL", year = 2015L)
)

abridged_list <- list()
truth_1y_list <- list()

for (t in targets) {
  base_dir <- file.path("data-raw/hmd", t$code)
  
  # Abridged Dx & Ex, collapsed to 85+ open interval
  dx <- read_hmd_wide(file.path(base_dir, "Deaths_5x1.txt"), "Dx")[Year == t$year]
  ex <- read_hmd_wide(file.path(base_dir, "Exposures_5x1.txt"), "Ex")[Year == t$year]
  ab <- merge(dx, ex, by = c("Year", "Age", "sex"))
  ab[, age_start := as.integer(sub("[-+].*$", "", Age))]
  ab[, age := pmin(age_start, 85L)]
  
  ab_85 <- ab[, .(Dx = sum(Dx, na.rm = TRUE), Ex = sum(Ex, na.rm = TRUE)),
              by = .(year = Year, sex = as.character(sex), age)]
  ab_85[, country := t$code]
  setcolorder(ab_85, c("country", "year", "sex", "age", "Dx", "Ex"))
  abridged_list[[length(abridged_list) + 1L]] <- ab_85
  
  # True 1x1 rates and life expectancy (ages 0 to 110+)
  mx1 <- read_hmd_wide(file.path(base_dir, "Mx_1x1.txt"), "mx_true")[Year == t$year]
  mx1[, age := as.integer(sub("\\+", "", Age))]
  
  lt_files <- c(Female = "fltper_1x1.txt", Male = "mltper_1x1.txt", Total = "bltper_1x1.txt")
  lt_1y <- rbindlist(lapply(names(lt_files), function(s) {
    d <- fread(file.path(base_dir, lt_files[[s]]), skip = 2, na.strings = ".")[Year == t$year]
    out <- d[, .(age = as.integer(sub("\\+", "", Age)), lx_true = lx / 1e5, ex_true = ex)]
    out[, sex := s]
    out
  }))
  
  truth <- merge(mx1[, .(sex = as.character(sex), age, mx_true)], lt_1y, by = c("sex", "age"))
  truth[, `:=`(country = t$code, year = t$year)]
  truth_1y_list[[length(truth_1y_list) + 1L]] <- truth
}

hmd_test_data <- list(
  std_input = std_input,
  abridged  = rbindlist(abridged_list),
  truth_1y  = rbindlist(truth_1y_list)
)

saveRDS(hmd_test_data, "inst/extdata/hmd_test_data.rds", compress = "xz")