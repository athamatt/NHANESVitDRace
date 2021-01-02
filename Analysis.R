# This file loads and analyzes NHANES data, 2011-14, usually using
# regression models, looking for relationships between serum 25(OH)D, race,
# socioeconomic variables, physiological and biochemical markers, and, as
# dependent variables, measures of cardiovascular risk, fracture risk, and
# periodontal health. See doi:xxxxx for details.

library(RNHANES) # Utility routines to load NHANES data
library(survey) # Handles stratified, weighted, multicluster samples
library(svyVGAM) # Count regression with survey data
library(xtable) # Generates Latex tables from R tables
library(Hmisc)  # Statistics functions
library(texreg) # Generates Latex tables from regression output
library(weights) # Histograms of data with weights
library(plyr) # rbind.fill() from here
library(MASS) # write.matrix() for writing out data
library(stringr) # String utility functions such str_detect

# Correct for single-non-certainty PSUs in a strata: conservative
# correction with the grand mean for strata mean
options(survey.lonely.psu="adjust")

# Utility functions
DEBUG_PRINT_ON <- TRUE;
tprint <- function(...) {if (DEBUG_PRINT_ON) print(...)}

TABLES_TO_PRINT <- c(1,2,3,13:46) # Which of the tables made here to output

## Variable recode functions. recode_<varname>

#From files DEMO_G,H:
# Set non-NHwhite ("3") to 0; non-NHBlack ("4") to 1; all else NA
recode_RIDRETH1 <- function(df) {ret <- dplyr::recode(df$RIDRETH1, "3" = 0, 
                                                      "4" = 1, .default = 2) 
                                 ret[ret > 1] <- NA
                                 return (as.integer(ret))}
#Sex: female ("2") to 0; male ("l") to 1.
recode_RIAGENDR <- function(df) dplyr::recode(df$RIAGENDR, "2" = 0, "1" = 1)
recode_RIDAGEYR <- function (df) return(df$RIDAGEYR) # age; tops at 80
recode_RIDEXMON <- function(df) return(df$RIDEXMON) # season, 1=winter; 2=summer
# Education: 1=less than 9th grade . . . 5 = college grad
recode_DMDEDUC2 <- function(df) {df$DMDEDUC2[(df$DMDEDUC2 == 7) | 
                                                (df$DMDEDUC2 == 9)] <- NA
                                 return(df$DMDEDUC2)}
recode_INDFMPIR <- function(df) df$INDFMPIR #fam inc/poverty-lvl. Capped at 5.

#From files SMQ_G,H (questionnaire)
#  SMQ040 == 3: nonsmoker
#  SMD641 = # of days smoked last month;
#  SMD650 = average # of cigs each round
recode_cig_smoking <- function(df) { # Return total cigs smoked last month
                       df$SMD641[(df$SMD641 == 77) |(df$SMD641 == 99)] <- NA
                       df$SMD650[(df$SMD650 == 777) | (df$SMD641 == 999)] <- NA
                       df$cig_smoking <- df$SMD641 * df$SMD650
                       df$cig_smoking[df$SMQ040 == 3] <- 0 # non-smoker
                       return(df$cig_smoking)
}

#From files TRIGLY_G,H (laboratory data)
# LDL cholesterol: LBDLDL no recoding (mg/dL): must use weight WTSAF2YR
recode_LBDLDL <- function(df) df$LBDLDL

#From files GLU_G,H (laboratory data)
## Fasting plasma glucose, no recoding: LBXGLU (mg/dL): must use weight WTSAF2YR
recode_LBXGLU <- function(df) df$LBXGLU

#From files BPX_G,H (examination data) 
#Blood pressure BPX_G,H: take means of 4 readings, all mm Hg
recode_blood_pressure_diastole <- function(df) rowMeans(cbind(df$BPXDI1, 
                                                              df$BPXDI2,
                                                          df$BPXDI3,df$BPXDI4),
                                                        na.rm = TRUE)
recode_blood_pressure_systole <- function(df) rowMeans(cbind(df$BPXSY1, 
                                                             df$BPXSY2,
                                                         df$BPXSY3,df$BPXSY4),
                                                       na.rm = TRUE)

#From files BMX_G,H (examination data)
#Body measures BMI, no recoding (kg/m^2)
recode_BMXBMI <- function(df) return(df$BMXBMI)

#From files DXXFRX_H (examination data): 
recode_hip_fracture_risk <- function(df) {
   df$hip_fracture_risk <- df$DXXFRAX3 # in DXXFRX: no previous hip fracture
   df$DXXPRVFX <- replace(df$DXXPRVFX, is.na(df$DXXPRVFX), 0) # Still NA
   # For previous fracture case use DXXFRAX1 as fracture score
   df$hip_fracture_risk[df$DXXPRVFX == 1] <- df$DXXFRAX1[df$DXXVFAST == 1]
   return(df$hip_fracture_risk)
}

# We use self-reported previous fracture to decide which risk variable to use
# for osteoporotic fractures as well.
recode_osteoporotic_fracture_risk <- function(df) {
   df$osteoporotic_fracture_risk <- df$DXXFRAX4
   df$DXXPRVFX <- replace(df$DXXPRVFX, is.na(df$DXXPRVFX), 0) # Required by R
   df$osteoporotic_fracture_risk[df$DXXPRVFX == 1] <- 
      df$DXXFRAX2[df$DXXPRVFX == 1]
   return(as.integer(round(df$osteoporotic_fracture_risk, 
                            digits = 1) * 10)) # Round to near 0.1 and scale
}

#From bone density files DXXSPN_H, DXXFEM_H (examination data)
# No recoding
recode_DXXOSBMD <- function(df) df$DXXOSBMD # Total spine BMD: only for > 40 yrs
recode_DXXOFBMD <- function(df) df$DXXOFBMD # Femur BMD: only for > 40 yrs

#File DXXAAC_H (examination data)
# No recoding: DXXAAC24
recode_DXXAAC24 <- function(df) df$DXXAAC24 # AAC24 score: only for > 40 yrs

#From files VID_G,H (laboratory data: WTMEC2YR weights)
# No recoding: LBXVIDMS (25 OH D2+D3 nmol/L): imputation used past det. limit
recode_LBXVIDMS <- function(df) df$LBXVIDMS

# From files: OHXPER_G,H (examintation data)

# Attachment loss (AL) from the many OHX..LA.. variables. Use means of the ALs
# from the various sites for the various teeth. Convert to integers by rounding
# and scaling.
recode_perio_att_loss <- function(df) {
   # The fields for attachment loss start with OHX, then have 2
   # digits, and then the substring LA. Search and extract these.
   tm <- df[str_detect(colnames(df), "OHX..LA")]
   tm[tm == 99] <- NA # 99 is the only invalid value
   rM <- rowMeans(tm, na.rm = TRUE) # Mean for all the teeth/sites
   rM <- rM - min(rM, na.rm = TRUE) # Adjust to min 0 for ease of handling
   return(as.integer(round(rM, digits = 1) * 10)) # Round to near 0.1 and scale
}


# Filter to obtain non-Hispanic Black/white data standalone. Expects variabe to
# be raw, that is, not recoded.
filter_in_nHBlack <- function(df) subset(df, recode_RIDRETH1(df) == 1)
filter_in_nHwhite <- function(df) subset(df, recode_RIDRETH1(df) == 0)

# Functions to print LR-output in Latex form using texreg
extract.svyglm <- function(model) {
   s <- summary(model)
   gof <- c(AIC(s$fit), logLik(s$fit), nobs(s$fit))
   gof.names <- c("AIC", "Log likelihood", "Num.\\ obs.")
   
   return(createTexreg(coef.names = rownames(s$coeftable),
                       coef = s$coeftable[, 1],
                       se = SE(model),
                       pvalues = s$coeftable[, 4],
                       gof.names = gof.names,
                       gof = gof
                       ))
}
setMethod("extract", signature = className("svy_vglm", "survey"),
           definition = extract.svyglm)
printLRoutput <- function(model,
                          caption,
                          coefnames,
                          dcolumn_on,
                          Fcount) {
   regobj <- extract(model)
   labelfull = paste0("Results:", Fcount$Tablen) # Table label

   # Now produce LATEX table code
   texreg(regobj, booktabs = TRUE, center = FALSE, dcolumn = dcolumn_on,
          custom.note = "***$\\!p < .001$; **$\\!p < .01$; *$\\!p < .05$",
          label = labelfull, caption = caption, caption.above = TRUE,
          longtable = TRUE, single.row = TRUE, custom.coef.names = coefnames,
          file = sprintf("%s%d.tex", "Table", Fcount$Tablen), digits = 3)
}

### The actual analysis code

# Merge side by side; for duplicated cols, retain name of one
mergetwo <- function (x, y) merge(x, y, by = "SEQN", suffixes = c("", ".y")) 

# Load a file across years. Join by attaching end to end.
load_multi_years <- function(file, years) {
   rbind.fill(nhanes_load_data(file, year = years, cache = TRUE))
}

# Plot the histogram of a weighted variable, with and without the 0 bin.
# Also, write the bin data (midpoint and count) as a 2-column .dat file
# var: variable name as a string.
# breaks: See help(hist). Set to number of bins to override default.
PlotHistDependentVar <- function(files_needed, years, var, weights = "WTMEC2YR",
                                 breaks = "Sturges")
{
   # Load all the files across the years into a list of frames
   dframes <- lapply(files_needed, load_multi_years, years)
   dff <- Reduce(mergetwo, dframes) # Merge into a single frame
   # Create a string containing the whole command dff$var <- recode_var(dff)
   exec_str <- paste0("dff$", var, " <- ", "recode_", var, "(dff)")
   tprint(exec_str)
   eval(parse(text = exec_str)) # Execute that command
   
   rdata <- dff[[var]]
   weightvals <- dff[[weights]]
   for (fname in c("raw", "no0s")) {
      # Produce a rough R draft figure of the data, with and without 0s
      hst_dat <- wtd.hist(rdata, col = "red", xlab = "", ylab = "", 
                          mgp = c(0,0,-0.5),tck = -0.01, ann = FALSE, 
                          weight = weightvals, breaks = breaks)
      title(paste0(var, " histogram ", fname), font = 2, line = 1)
      mtext(text = paste0(var, " raw score"), side = 1, line = 1.2, font = 2)
      mtext(text = "Count", side = 2, line = 1.2, font = 2)
      # Dump the bin data to a .dat file for pgfplot; gets us a finished pic
      write.matrix(matrix(c(hst_dat$breaks, hst_dat$counts), 
                          nrow = length(hst_dat$breaks), ncol = 2, 
                          byrow = FALSE), 
                   file = paste0(var, fname, ".dat"), sep = " ")
      selcols <- (rdata > 0)
      rdata<- rdata[selcols] # Zero-removed data next pass
      weightvals <- weightvals[selcols] # Adjust weight vector to match
   }
}

# Figure out the sample weight to use. We default to WTMEC2YR, the examined
# samples, as our dvars are vars like AAC24, periodontal stuff, DEXA stuff, and
# so on. In some cases we have fasting-based measures such as plasma glucose
# and LDL/triglycerides. For that still another subsample weight, WTSAF2YR, is
# needed. The narrowest measure must be used: WTSAF2YR over WTMEC2YR over 
# WTINT2YR.
get_weights <- function(vars) {
   weights <- "WTMEC2YR"
   fasting_var <- "LBXGLU|LBDGLUSI|LBXIN|LBDINSI|LBXTR|LBDTRSI|LBDLDL|LBDLDLSI";
   # Convert ahove string to "^LBXGLU$|^LBDGLUSI$|....|^LBDLDLSI$"
   # That ensures we match the whole variable name and not a substring
   fasting_var <- gsub("|", "$|^", fasting_var, fixed = TRUE)
   fasting_var <- gsub("$", "$", gsub("^", "^", fasting_var))
   if (any(grepl(fasting_var, vars))) {
      weights <- "WTSAF2YR"
      tprint("Using WTSAF2YR")
   }
   return (weights)
}

# Run one analysis. Specify the files we need to load, the years of interest,
# the dependent variables as a string list, the independent variable, and, 
# optionally, whether to use GLM or VGLM, what count regression to use, and
# any filtering to perform before the recoding and the analysis.
run_analysis <- function(files_needed, years, ivar_strings, dvar, 
                         interactions = NULL, use_glm = FALSE, 
                         vglm_fn = zinegbinomialff, filter_fn = NULL) {
   # Load all the files across the years into a list of frames
   dframes <- lapply(files_needed, load_multi_years, years)
   dff <- Reduce(mergetwo, dframes) # Merge into a single frame

   if (!is.null(filter_fn)) { # Invoke data filter if argument passed in
      tprint("Calling filter function: ")
      dff <- filter_fn(dff)
   }
   
   # Recode the independent variables
   for (i in ivar_strings) {
      tprint(i)
      # Create a string containing the whole command dff$var <- recode_var(dff)
      exec_str <- paste0("dff$", i, " <- ", "recode_", i, "(dff)")
      tprint(exec_str)
      eval(parse(text = exec_str)) # Execute that command
   }
   # Now recode the DV
   exec_str <- paste0("dff$", dvar, " <- ", "recode_", dvar, "(dff)")
   tprint(exec_str)
   eval(parse(text = exec_str))

   # Find the appropriate weights. Then adjust the weights based on the number 
   # of cycles, that is, the length of the years vector.
   weights <- get_weights(c(ivar_strings, dvar))
   dff[[weights]] <- dff[[weights]] / length(years)
   
   # Make the survey design object with the RNHANES functions.
   # Strata, VPSU and weights specified.
   NHANES_des <- svydesign(data = dff, ids = ~SDMVPSU, strata = ~SDMVSTRA, 
                           weights = as.formula(paste0("~", weights)), 
                           nest = TRUE)
   
   # Create the GLM. The model formula is dynamically created as a string
   ivar_strings <- c(ivar_strings, interactions)
   var_list<-paste0(unlist(ivar_strings), collapse = "+")
   var_list<-paste0(dvar, "~",  var_list)
   if (use_glm) {
      tprint("GLM: ")
      tprint(var_list)
      reg_result <- svyglm(as.formula(var_list), design = NHANES_des)
   } else {
      tprint("VGLM: ")
      tprint(var_list)
      reg_result <- svy_vglm(formula = as.formula(var_list), family = vglm_fn(), 
                             design = NHANES_des, crit = "coef", maxit = 200)
   }
   return(reg_result)
}

summarizePrint <- function(rslt, title, Fcount, 
                           tables_to_print = TABLES_TO_PRINT)
{
   if (Fcount$LRn %in% tables_to_print) {
      printLRoutput(rslt, title, coefnames = names(rslt$coef), 
                    dcolumn_on = TRUE, Fcount)
      Fcount$Tablen <- Fcount$Tablen + 1
   }
   Fcount$LRn <- Fcount$LRn + 1
   return(Fcount)
}

Fcount <- data.frame(LRn = 1, Tablen = 1) #LRn regression, Tablen table, counts

#################### Cardio vascular risk analysis #############################

# Histograms of dependent variable
PlotHistDependentVar(c("DEMO", "DXXAAC"), "2013-2014", "DXXAAC24")

# Start the analyses
# AAC vs race
rslt <- run_analysis(files_needed = c("DXXAAC", "DEMO"), 
                     years <- c("2013-2014"), ivar_strings <- c("RIDRETH1"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race", Fcount)

# AAC vs vit-D (25 hydroxy)
rslt <- run_analysis(files_needed = c("DXXAAC", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs vit D", Fcount)


# AAC vs race and vit-D (25 hydroxy) and interaction
rslt <- run_analysis(files_needed = c("DXXAAC", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("DXXAAC24"),
                     interactions <- c("LBXVIDMS*RIDRETH1"))
Fcount <- summarizePrint(rslt, "AAC vs vit D, race, interaction", Fcount)

# AAC vs vit D, edu. level, and interaction
rslt <- run_analysis(files_needed = c("DXXAAC", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "DMDEDUC2"), 
                     dvar <- c("DXXAAC24"),
                     interactions <- c("LBXVIDMS*DMDEDUC2"))
Fcount <- summarizePrint(rslt, "AAC vs vit D, edu., interaction", Fcount)


# AAC vs race, LDL
rslt <- run_analysis(files_needed = c("DXXAAC", "TRIGLY", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings = c("RIDRETH1", "LBDLDL"), 
                     dvar = c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, LDL", Fcount)

# AAC vs race, diastolic BP
rslt <- run_analysis(files_needed = c("DXXAAC", "BPX", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_diastole"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, diastolic BP", Fcount)

# AAC vs race, systolic BP
rslt <- run_analysis(files_needed = c("DXXAAC", "BPX", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, systolic BP", Fcount)

# AAC vs race, fasting plasma glucose
rslt <- run_analysis(files_needed = c("DXXAAC", "GLU", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "LBXGLU"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, blood glucose", Fcount)

# Combination: AAC vs race plus health indicators: LDL, systolic BP, vit D, 
# edu level
rslt <- run_analysis(files_needed = c("DXXAAC", "TRIGLY","BPX", "DEMO", "VID"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole", 
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, sys. BP, LDL, edu., vit D", 
                          Fcount)

# Combination: AAC vs race plus BMI
rslt <- run_analysis(files_needed = c("DXXAAC", "DEMO", "VID", "BMX"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "BMXBMI", "LBXVIDMS", 
                                       "RIDEXMON"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs race, BMI", Fcount)

# Combination: AAC vs vit D, race plus cigarette smoking
rslt <- run_analysis(files_needed = c("DXXAAC", "DEMO", "VID", "SMQ"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "cig_smoking", "LBXVIDMS", 
                                       "RIDEXMON"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs vit D, race, smoking", Fcount)

# Combination: AAC vs race plus health indicators: LDL, systolic BP, vit D, 
# edu level, BMI
rslt <- run_analysis(files_needed = c("DXXAAC", "TRIGLY","BPX", "DEMO", "VID", 
                                      "BMX"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole", 
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, "AAC vs vit D, LDL, edu., race, sys. BP, BMI", 
                          Fcount)

# AAC vs race, vit D, smoking, income, edu., BMI, LDL, glucose, sys. BP, age, sex
rslt <- run_analysis(files_needed = c("DXXAAC", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("DXXAAC24"))
Fcount <- summarizePrint(rslt, paste0("AAC vs race, vit D, smoking, income, ", 
                                       "edu., BMI, LDL, glucose, sys. BP, ",
                                       "age, sex"),
                          Fcount)

# AAC vs vit-D (25 hydroxy): non-Hispanic Black only
rslt <- run_analysis(files_needed = c("DXXAAC", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXAAC24"), filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt, "AAC vs vit D: nHBlacks only", Fcount)

# AAC vs vit-D (25 hydroxy): non-Hispanic white only
rslt <- run_analysis(files_needed = c("DXXAAC", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXAAC24"), filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt, "AAC vs vit D: nHwhites only", Fcount)

######################## Bone density/fracture risk analyses ###################

PlotHistDependentVar(c("DEMO", "DXXSPN"), "2013-2014", "DXXOSBMD")

# spine density vs vit D (25 OH D)
rslt <- run_analysis(files_needed = c("DXXSPN", "VID", "DEMO"), 
                     years <- c("2013-2014"),
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOSBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, "SBMD vs vit D", Fcount)

# spine density vs vit D vs race
rslt <- run_analysis(files_needed = c("DXXSPN", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("DXXOSBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, "SBMD vs vit D, race", Fcount)

# spine density vs vit D/race/interaction
rslt <- run_analysis(files_needed = c("DXXSPN", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("DXXOSBMD"), 
                     use_glm = TRUE, interactions = c("LBXVIDMS*RIDRETH1"))
Fcount <- summarizePrint(rslt, "SBMD vs vit D, race, interaction", Fcount)

# spine density vs race, vit D, smoking, income, edu., BMI, LDL, glucose, 
# sys. BP, age, sex
rslt <- run_analysis(files_needed = c("DXXSPN", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("DXXOSBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, paste0("SBMD vs race, vit D, smoking, income, ", 
                                       "edu., BMI, LDL, glucose, sys. BP, ",
                                       "age, sex"),
                          Fcount)

# spine density vs vit D: nHBlacks only
rslt <- run_analysis(files_needed = c("DXXSPN", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOSBMD"), use_glm = TRUE,
                     filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt, "SBMD vs vit D: nHBlacks only", Fcount)

# spine density vs vit D: nHwhites only
rslt <- run_analysis(files_needed = c("DXXSPN", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOSBMD"), use_glm = TRUE,
                     filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt, "SBMD vs vit D: nHwhites only", Fcount)

PlotHistDependentVar(c("DEMO", "DXXFEM"), "2013-2014", "DXXOFBMD")

# femur density vs vit D (25 OH D)
rslt <- run_analysis(files_needed = c("DXXFEM", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, "FBMD vs vit D", Fcount)

# femur density vs vit D vs race
rslt <- run_analysis(files_needed = c("DXXFEM", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, "FBMD vs vit D, race", Fcount)

# femur density vs vit D/race/interaction
rslt <- run_analysis(files_needed = c("DXXFEM", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE, 
                     interactions = c("LBXVIDMS*RIDRETH1"))
Fcount <- summarizePrint(rslt, "FBMD vs vit D, race, interaction", Fcount)

# femur density vs race, vit D, smoking, income, edu., BMI, LDL, glucose, 
# sys. BP, age, sex
rslt <- run_analysis(files_needed = c("DXXFEM", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE)
Fcount <- summarizePrint(rslt, paste0("FBMD vs race, vit D, smoking, income, ", 
                                       "edu., BMI, LDL, glucose, sys. BP, ",
                                       "age, sex"),
                          Fcount)

# femur density vs vit D: nHBlacks only
rslt <- run_analysis(files_needed = c("DXXFEM", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE,
                     filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt, "FBMD vs vit D: nHBlacks only", Fcount)

# femur density vs vit D: nHwhites only
rslt <- run_analysis(files_needed = c("DXXFEM", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("DXXOFBMD"), use_glm = TRUE,
                     filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt, "FBMD vs vit D: nHwhites only", Fcount)


PlotHistDependentVar(c("DEMO", "DXXFRX", "DXXVFA"), "2013-2014",
                     "hip_fracture_risk")

# Hip fracture risk vs Vit D
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("hip_fracture_risk"), 
                     use_glm = TRUE)
Fcount <- summarizePrint(rslt, "Hip fracture risk vs vit D", Fcount)

# Hip fracture risk vs Vit D, race, and interaction
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("hip_fracture_risk"), 
                     use_glm = TRUE, interactions = c("LBXVIDMS*RIDRETH1"))
Fcount <- summarizePrint(rslt, "Hip fracture risk vs vit D, race, interaction", 
                          Fcount)

# Hip fracture risk vs vit D, race, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "DEMO", "BPX", "VID", 
                                      "TRIGLY", "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("hip_fracture_risk"), 
                     use_glm = TRUE)
Fcount <- summarizePrint(rslt, paste0("Hip fracture risk vs vit D, race, ",
                                       "systolic BP, LDL, edu., BMI, ",
                                       "smoking, glucose, income ratio, ",
                                       "age, sex "), 
                          Fcount)

# Hip fracture risk vs Vit D: nHBlacks only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("hip_fracture_risk"), use_glm = TRUE, 
                     filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt, "Hip fracture risk vs vit D: nHBlacks only", 
                          Fcount)

# Hip fracture risk vs Vit D: nHwhites only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("hip_fracture_risk"), use_glm = TRUE, 
                     filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt, "Hip fracture risk vs vit D: nHwhites only", 
                          Fcount)


PlotHistDependentVar(c("DEMO", "DXXFRX", "DXXVFA"), "2013-2014",
                     "osteoporotic_fracture_risk")

# Osteoporotic fracture risk vs Vit D
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt, "Osteo. fracture risk vs vit D", Fcount)

# Osteoporotic fracture risk vs Vit D, race, and interaction
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON", "RIDRETH1"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial, 
                     interactions = c("LBXVIDMS*RIDRETH1"))
Fcount <- summarizePrint(rslt, paste0("Osteo. fracture risk vs vit D, race, ", 
                                       "interaction"), Fcount)

# Osteo fracture risk vs vit D, race, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "DEMO", "BPX", "VID", 
                                      "TRIGLY", "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt,
                          paste0("Osteo. frac. risk vs vit D, race, ",
                          "syst. BP, LDL, edu., BMI, smoking, glucose, ",
                          "income ratio, age, sex"), Fcount)

# Osteoporotic fracture risk vs Vit D: nHBlacks only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial, filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt, "Osteo. fracture risk vs vit D: nHBlacks only", 
                          Fcount)

# Osteoporotic fracture risk vs Vit D: nHwhites only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "VID", "DEMO"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial, filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt, "Osteo. fracture risk vs vit D: nHwhites only", 
                          Fcount)

# Osteo fracture risk vs vit D, race, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex: nHBlacks only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "DEMO", "BPX", "VID", 
                                      "TRIGLY", "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial, filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt,
                          paste0("Osteo. frac. risk vs vit D, syst. BP, LDL, ",
                                 "edu., BMI, smoking, glucose, ",
                                 "income ratio, age, sex: nHBlack"), Fcount)

# Osteo fracture risk vs vit D, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex: nHwhites only
rslt <- run_analysis(files_needed = c("DXXFRX", "DXXVFA", "DEMO", "BPX", "VID", 
                                      "TRIGLY", "BMX", "SMQ", "GLU"), 
                     years <- c("2013-2014"), 
                     ivar_strings <- c("blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("osteoporotic_fracture_risk"), 
                     vglm_fn = negbinomial, filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt,
                          paste0("Osteo. frac. risk vs vit D, syst. BP, LDL, ",
                                 "edu., BMI, smoking, glucose, ",
                                 "income ratio, age, sex: nHwhite"), Fcount)

########################## Dental health and vit-D #############################

PlotHistDependentVar(c("DEMO", "OHXPER"), c("2011-2012", "2013-2014"), 
                     "perio_att_loss")

# Periodontitis vs Vit D
rslt <- run_analysis(files_needed = c("OHXPER", "VID", "DEMO"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("perio_att_loss"), 
                     use_glm = FALSE, vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt, "AL vs vit D", Fcount)

# Periodontitis vs race
rslt <- run_analysis(files_needed = c("OHXPER", "DEMO"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("RIDRETH1"), dvar <- c("perio_att_loss"), 
                     use_glm = FALSE, vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt, "AL vs race", Fcount)

# Periodontitis vs vit D and race
rslt <- run_analysis(files_needed = c("OHXPER", "VID", "DEMO"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("perio_att_loss"), 
                     use_glm = FALSE, vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt, "AL vs vit D, race", Fcount)

# Periodontitis vs vit D, race, and interaction
rslt <- run_analysis(files_needed = c("OHXPER", "VID", "DEMO"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "LBXVIDMS", "RIDEXMON"), 
                     dvar <- c("perio_att_loss"), 
                     interactions = c("RIDRETH1*LBXVIDMS"),
                     use_glm = FALSE, vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt, "AL vs vit D, race, interaction", Fcount)

# Periodontitis vs vit D, race, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex
rslt <- run_analysis(files_needed = c("OHXPER", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("RIDRETH1", "blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("perio_att_loss"), vglm_fn = negbinomial)
Fcount <- summarizePrint(rslt,
                          paste0("AL vs vit D, race, systolic BP, LDL, ",
                          "edu., BMI, smoking, glucose, income ratio, ",
                          "age, sex"), 
                          Fcount)

# Periodontitis vs vit D, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex: nHBlack only
rslt <- run_analysis(files_needed = c("OHXPER", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("perio_att_loss"), vglm_fn = negbinomial,
                     filter_fn = filter_in_nHBlack)
Fcount <- summarizePrint(rslt,
                          paste0("AL vs vit D, systolic BP, LDL, ",
                                 "edu., BMI, smoking, glucose, ",
                                 " income ratio, age, sex: nHBlack only"), 
                          Fcount)

# Periodontitis vs vit D, systolic BP, LDL, edu., BMI, smoking, 
# glucose, income ratio, age, sex; nHwhite only
rslt <- run_analysis(files_needed = c("OHXPER", "DEMO", "BPX", "VID", "TRIGLY",
                                      "BMX", "SMQ", "GLU"), 
                     years <- c("2011-2012", "2013-2014"), 
                     ivar_strings <- c("blood_pressure_systole",
                                       "LBDLDL", "DMDEDUC2", "LBXVIDMS", 
                                       "RIDEXMON", "BMXBMI", "cig_smoking",
                                       "LBXGLU", "INDFMPIR", "RIDAGEYR", 
                                       "RIAGENDR"), 
                     dvar <- c("perio_att_loss"), vglm_fn = negbinomial,
                     filter_fn = filter_in_nHwhite)
Fcount <- summarizePrint(rslt,
                          paste0("AL vs vit D, systolic BP, LDL, ",
                                 "edu., BMI, smoking, glucose, ",
                                 " income ratio, age, sex: nHwhite only"), 
                          Fcount)
