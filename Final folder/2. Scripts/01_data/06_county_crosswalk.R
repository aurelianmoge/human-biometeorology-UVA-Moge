# =============================================================================
# 06_county_crosswalk.R - county / independent city -> climate region, and regional
# population denominators.
#
# The 12 workbooks partition Virginia by the decedent's county of residence
# (Davis et al. 2023, s.2.1.3; DeGuzman et al. 2026, Fig. 1), but the files do not list
# the counties.  Rules, in order:
#   1. Davis et al. (2025, Int J Biometeorol, Suppl. Table 2) roster for the 7 metro
#      regions (IAD CHO LYH PHF ORF RIC ROA); its death totals match our workbooks
#      within 0.3% (ROA exactly);
#   2. localities enclosed by / enclosing roster localities but omitted from the printed
#      list: Roanoke County, Covington, Radford -> ROA; Manassas, Manassas Park -> IAD;
#   3. Davis et al. (2023): the former Martinsville region (Southside) was merged into
#      Emporia (EMV): Henry, Patrick, Martinsville, Pittsylvania, Danville, Franklin;
#   4. every other locality -> nearest of the 5 non-metro stations (EMV EZF OKV SHD VJI).
# Rules 3-4 are inferred and flagged; to be confirmed with the authors (Fig. 1 map).
# Check: indirectly age-standardised mortality ratio (SMR) of each workbook region.
# Input : 1. Data/reference/metadata/county_region_crosswalk.csv, stations.csv;
#         1. Data/derived/county_joint_populations.csv, regional_panel.csv
# Output: 1. Data/derived/county_region_crosswalk.csv, regional_populations.csv;
#         3. Outputs/tables/01_data/crosswalk_validation.csv
# =============================================================================
source("2. Scripts/00_setup/config.R")
STAGE <- "01_data"
cw  <- fread(file.path(META, "county_region_crosswalk.csv")); stn <- fread(file.path(META, "stations.csv"))
pop <- fread(file.path(DERIVED, "county_joint_populations.csv"))
cw <- cw[FIPS %in% unique(pop$FIPS)]                    # Bedford city merged into Bedford County in the Census file

hav <- function(la1, lo1, la2, lo2) { r <- pi / 180; a <- sin((la2 - la1) * r / 2)^2 + cos(la1 * r) * cos(la2 * r) * sin((lo2 - lo1) * r / 2)^2; 12742 * asin(sqrt(a)) }
rest <- stn[Region %in% c("EMV", "EZF", "OKV", "SHD", "VJI")]
enclave   <- c("Roanoke County" = "ROA", "Covington City" = "ROA", "Radford City" = "ROA", "Manassas City" = "IAD", "Manassas Park City" = "IAD")
southside <- c("Henry County", "Patrick County", "Martinsville City", "Pittsylvania County", "Danville City", "Franklin County")
cw[, nearest_nonmetro := rest$Region[which.min(hav(Latitude, Longitude, rest$Latitude, rest$Longitude))], by = FIPS]
cw[, `:=`(Region = nearest_nonmetro, Basis = "inferred: nearest non-metro station")]
cw[Locality %in% southside, `:=`(Region = "EMV", Basis = "inferred: Davis 2023 (Martinsville region merged into EMV)")]
cw[Locality %in% names(enclave), `:=`(Region = enclave[Locality], Basis = "enclave of a Davis 2025 roster locality")]
cw[Reference_region != "", `:=`(Region = Reference_region, Basis = "Davis 2025 Suppl. Table 2")]
fwrite(cw[, .(FIPS, Locality, Region, Basis, Nearest_any_station = Nearest_sensitivity)], file.path(DERIVED, "county_region_crosswalk.csv"))

# ---- regional populations by group and year ----------------------------------------------------
pg <- rbind(pop[, .(Group = "All", P = sum(Population)), by = .(Year, FIPS)], pop[, .(P = sum(Population)), by = .(Year, FIPS, Group = Sex)],
            pop[Race %in% c("Black", "White"), .(P = sum(Population)), by = .(Year, FIPS, Group = Race)], pop[, .(P = sum(Population)), by = .(Year, FIPS, Group = Age)])
rp <- merge(pg, cw[, .(FIPS, Region)], by = "FIPS")[, .(Population = sum(P)), by = .(Region, Year, Group)]
fwrite(rp, file.path(DERIVED, "regional_populations.csv"))

# ---- validation: indirectly age-standardised mortality ratio per workbook region ---------------------
panel <- fread(file.path(DERIVED, "regional_panel.csv"))
d <- melt(panel[, c("Region", "Date", GROUPS), with = FALSE], id.vars = c("Region", "Date"), variable.name = "Group", value.name = "deaths")[
  , Year := year(as.IDate(Date))][, .(deaths = sum(deaths, na.rm = TRUE) * 365.25 / .N), by = .(Region, Year, Group)]   # annualised (missing days)
ages <- c("Under50", "Age50_74", "Age75plus")
pa <- pop[, .(P = sum(Population)), by = .(Year, FIPS, Group = Age)]
srate <- merge(d[Group %in% ages, .(D = sum(deaths)), by = .(Year, Group)], pa[, .(P = sum(P)), by = .(Year, Group)])[, .(Year, Group, rate = D / P)]
smr <- function(map) { rpa <- merge(pa, map, by = "FIPS")[, .(P = sum(P)), by = .(Region, Year, Group)]
  ex <- merge(rpa, srate, by = c("Year", "Group"))[, .(expected = sum(P * rate)), by = Region]
  merge(d[Group == "All", .(observed = sum(deaths)), by = Region], ex, by = "Region")[, SMR := observed / expected] }
v <- merge(smr(cw[, .(FIPS, Region)]), smr(cw[, .(FIPS, Region = Nearest_sensitivity)])[, .(Region, SMR_nearest_any = SMR)], by = "Region")
v <- merge(v, cw[, .(localities = .N, inferred = sum(grepl("inferred", Basis))), by = Region], by = "Region")
write_tab(v[order(SMR)], STAGE, "crosswalk_validation")
cat("SMR range", paste(round(range(v$SMR), 2), collapse = "-"), "| SD log SMR", round(sd(log(v$SMR)), 3),
    "(nearest-station crosswalk:", round(sd(log(v$SMR_nearest_any)), 3), ")\n")
