

preprocess_xfplate <- function(xf){


  # Use our xf list from read_xfplate() with all the necessary Seahorse data to fill this data tibble.

  xf_raw_pr <- preprocess_xf_raw(xf$raw,
                                 xf$pHcal,
                                 xf$inj,
                                 xf$assayinfo,
                                 xf$buffer,
                                 xf$norm,
                                 xf$flagged)

  xf_rate_pr <- preprocess_xf_rate(xf$rate,
                                   xf$norm,
                                   xf$flagged)

  xf_plate_pr <- xf_raw_pr %>%
    dplyr::group_by(plate_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(filepath_seahorse = list(tibble::tibble(
                  directory_path = dirname(as.character(xf$filepath_seahorse)),
                  base_name = basename(as.character(xf$filepath_seahorse)),
                  full_path = xf$filepath_seahorse
                )),
                  date = xf$assayinfo$date_run,
                  assay_info = list(tibble::tibble(xf$assayinfo)),
                  rate_data = list(tibble::tibble(xf_rate_pr)),
                  injection_info = list(tibble::tibble(xf$inj))) %>%
    dplyr::select(plate_id, filepath_seahorse, date, assay_info, injection_info,
                  raw_data = data, rate_data)


  return(xf_plate_pr)
}

preprocess_xf_raw <- function(xf_raw,
                              xf_pHcal,
                              xf_inj,
                              xf_assayinfo,
                              xf_buffer,
                              xf_norm,
                              xf_flagged) {


  # convert the original integer column to integers again, instead of double
  xf_raw_pr <- xf_raw %>%
    tibble::as_tibble() %>%
    dplyr::mutate(across(c(Measurement,
                           Tick,
                           `O2 Light Emission`,
                           `O2 Dark Emission`,
                           `O2 Ref Light`,
                           `O2 Ref Dark`,
                           `pH Light`,
                           `pH Dark`,
                           `pH Ref Light`,
                           `pH Ref Dark`),
                         as.integer))

  # rename columns
  xf_raw_pr <- rename_columns(xf_raw_pr)

  # convert time column
  xf_raw_pr <- convert_timestamp(xf_raw_pr)

  # correct pH_em_corr
  xf_raw_pr$pH_em_corr_corr <- correct_pH_em_corr(xf_raw_pr$pH_em_corr,
                                                  xf_pHcal$pH_cal_em,
                                                  xf_assayinfo$pH_targetEmission[1])


  # calculate backgrounds and join
  background <- calc_background(xf_raw_pr)

  xf_raw_pr <- xf_raw_pr %>%
    dplyr::left_join(background, by = c("measurement"), relationship = "many-to-many")

  # add injection info
  xf_raw_pr <- dplyr::left_join(xf_raw_pr, xf_inj, by = "measurement")

  #add plate_id to df
  xf_raw_pr$plate_id <- xf_assayinfo$plate_id

  #add norm_info
  xf_raw_pr <- xf_raw_pr %>% dplyr::left_join(xf_norm, by = c("well"))

  #add bufferfactor
  xf_raw_pr <- xf_raw_pr %>% dplyr::left_join(xf_buffer, by = c("well"))

  #add flag well columnn
  xf_raw_pr$flagged_well <- FALSE
  xf_raw_pr$flagged_well[xf_raw_pr$well %in% xf_flagged] <- TRUE

  # select columns that are needed
  xf_raw_pr <- xf_raw_pr %>% dplyr::select(
    plate_id, well, measurement, tick, timescale, minutes, group, interval, injection,
    O2_em_corr, pH_em_corr, O2_mmHg, pH, pH_em_corr_corr, O2_em_corr_bkg,
    pH_em_corr_bkg, O2_mmHg_bkg, pH_bkgd, pH_em_corr_corr_bkg, bufferfactor, cell_n, flagged_well
  )

  return(xf_raw_pr)

  }

preprocess_xf_rate <- function(xf_rate,
                               xf_norm,
                               xf_flagged){
  #add norm_info to rate data
  OCR_from_excel <- xf_rate %>% dplyr::left_join(xf_norm, by = c("well"))

  OCR_from_excel$flagged_well <- FALSE
  OCR_from_excel$flagged_well[OCR_from_excel$well %in% xf_flagged] <- TRUE

  return(OCR_from_excel)
}


rename_columns <- function(xf_raw_pr) {

  # change column names into terms without spaces
  colnames(xf_raw_pr) <- c(
    "measurement", "tick", "well", "group",
    "time", "temp_well", "temp_env", "O2_isvalid", "O2_mmHg",
    "O2_light", "O2_dark", "O2ref_light", "O2ref_dark",
    "O2_em_corr", "pH_isvalid", "pH", "pH_light", "pH_dark",
    "pHref_light",
    "pHref_dark", "pH_em_corr", "interval"
  )

  return(xf_raw_pr)
}

convert_timestamp <- function(xf_raw_pr) {

  # first make sure that the data is sorted correctly
  xf_raw_pr <- dplyr::arrange(xf_raw_pr, tick, well)

  # add three columns to df (totalMinutes, minutes and time) by converting the timestamp into seconds
  xf_raw_pr$time <- as.character((xf_raw_pr$time))
  times <- strsplit(xf_raw_pr$time, ":")
  xf_raw_pr$totalMinutes <- sapply(times, function(x) {
    x <- as.numeric(x)
    x[1] * 60 + x[2] + x[3] / 60
  })
  xf_raw_pr$minutes <- xf_raw_pr$totalMinutes - xf_raw_pr$totalMinutes[1] # first row needs to be first timepoint!
  xf_raw_pr$timescale <- round(xf_raw_pr$minutes * 60)

  return(xf_raw_pr)
}

correct_pH_em_corr <- function(pH_em_corr, pH_cal_em, pH_targetEmission){

  correct_pH_em_corr <- (pH_targetEmission / pH_cal_em) * pH_em_corr

}


calc_background <- function(xf_raw_pr){

  background <- xf_raw_pr %>%
    dplyr::select(group, well, measurement, timescale, O2_em_corr,
           pH_em_corr, O2_mmHg, pH, pH_em_corr_corr) %>%
    dplyr::filter(group == "Background") %>%
    dplyr::reframe(
      measurement,
      O2_em_corr_bkg = mean(O2_em_corr),
      pH_em_corr_bkg = mean(pH_em_corr),
      O2_mmHg_bkg = mean(O2_mmHg),
      pH_bkgd = mean(pH),
      pH_em_corr_corr_bkg = mean(pH_em_corr_corr)
    )

  return(background)
}

