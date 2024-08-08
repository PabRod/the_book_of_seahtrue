
read_xfplate <- function(filepath_seahorse) {

    # read data
    xf_raw <- get_xf_raw(filepath_seahorse)
    xf_rate <- get_xf_rate(filepath_seahorse) #outputs list of 2
    xf_norm <- get_xf_norm(filepath_seahorse) #outputs list of 2
    xf_buffer <- get_xf_buffer(filepath_seahorse)
    xf_inj <- get_xf_inj(filepath_seahorse)
    xf_pHcal <- get_xf_pHcal(filepath_seahorse)
    xf_O2cal <- get_xf_O2cal(filepath_seahorse)
    #xf_flagged <- get_xf_flagged(filepath_seahorse)
    xf_assayinfo <- get_xf_assayinfo(filepath_seahorse,
                                     norm_available = xf_norm[[2]],
                                     xls_ocr_backgroundcorrected = xf_rate[[2]])
    xf_norm <- xf_norm[[1]]
    xf_rate <- xf_rate[[1]]

    # make the output list
    xf <- list(
      raw = xf_raw,
      rate = xf_rate,
      assayinfo = xf_assayinfo,
      inj = xf_inj,
      pHcal = xf_pHcal,
      O2cal = xf_O2cal,
      norm = xf_norm,
      flagged = "empty",
      buffer = xf_buffer,
      filepath_seahorse = filepath_seahorse
    )

    logger::log_info(glue::glue("Parsing all collected seahorse information from file: {filepath_seahorse}"))

    return(xf)

}

get_xf_raw <- function(filepath_seahorse){

   try(
    xf_raw <- readxl::read_excel(filepath_seahorse,
                         sheet = "Raw",
                         col_types = c("numeric", # Measurment
                                       "numeric", # Tick
                                       "text", # Well
                                       "text", # Group
                                       "text", # TimeStamp
                                       "numeric", # Well Temperature
                                       "numeric", # Environment Temperature
                                       "text", # O2 is Valid
                                       "numeric", # O2 (mmHg)
                                       "numeric", # O2 Light Emission
                                       "numeric", # O2 Dark Emission
                                       "numeric", # O2 Ref Light
                                       "numeric", # O2 Ref Dark
                                       "numeric", # O2 Corrected Em.
                                       "text", # pH Is Valid
                                       "numeric", # pH
                                       "numeric", # pH Light
                                       "numeric", # pH Dark
                                       "numeric", # pH Ref Light
                                       "numeric",# pH Ref Dark
                                       "numeric" # pH Corrected Em.
                         ))
   )

   if (exists("xf_raw")) {
     cli::cli_alert_success("Succesfully read Raw information from data sheet.")
   } else {
     cli::cli_abort("An error occured during the data collection from 'Raw'sheet.")
   }


    logger::log_info("Finished collecting data from 'Raw' sheet.")

    return(xf_raw)

      }

get_xf_norm <- function(filepath_seahorse){

    norm_info <- get_platelayout_data(filepath_seahorse,
                                      my_sheet = "Assay Configuration",
                                      my_range = "B84:N92",
                                      my_param = "cell_n")


    if (sum(is.na(norm_info$cell_n)) >90){
      norm_available <- FALSE
    } else {
      norm_available <- TRUE}

    xf_norm <- list(norm_info, norm_available)


    if (exists("xf_norm")) {
      cli::cli_alert_success("Succesfully read normalisation info from 'Assay configuration' from data sheet.")
    } else {
      stop("An error occured during the data collection from 'Assay Configuration' sheet.")
    }


    logger::log_info("Finished collecting normalisation info from 'Assay configuration' sheet.")

    return(xf_norm)
}

get_originalRateTable <- function(filepath_seahorse){

  logger::log_info("Collecting OCR from 'Rate' sheet.")

  original_rate_df <- readxl::read_excel(filepath_seahorse, sheet = "Rate")

  # because rate data can be either background corrected or not this should be checked first
  # first verify whether a "Background" group exists in the  original_rate_df


  if ("Background" %in% {original_rate_df$Group %>% unique()}) {

    logger::log_info("A background group was found in the RATE sheet")


    check_background <- original_rate_df %>%
      dplyr::filter(Group == "Background") %>%
      dplyr::select(OCR) %>%
      dplyr::reframe(mean = mean(OCR)) %>%
      dplyr::pull(mean)

    if (check_background == 0) {
      corrected_allready <- TRUE
    } else {
      corrected_allready <-  FALSE
    }

  } else {

    #in case when there is no Background group we work with the original data
    # that is in the input file "Rate" sheet
    # please note that there will be warning logged, but the columns will be
    # labeled incorrectly as if the data is background corrected

    logger::log_info("WARNING: no background group was found in the 'Rate' sheet")

    corrected_allready <-  TRUE

  }

  if (corrected_allready == TRUE){
    colnames(original_rate_df) <-
      c("measurement","well", "group",
        "time_wave", "OCR_wave_bc",
        "ECAR_wave_bc", "PER_wave_bc")
    original_rate_df <- original_rate_df %>%
      dplyr::mutate(OCR_wave = 0, ECAR_wave = 0)

    original_rate_df <- original_rate_df %>%
      dplyr::select(measurement, well, group,
                    time_wave, OCR_wave, OCR_wave_bc,
                    ECAR_wave, ECAR_wave_bc)

  } else{
    colnames(original_rate_df) <-
      c("measurement","well", "group",
        "time_wave", "OCR_wave",
        "ECAR_wave", "PER_wave")

    #do background substraction forr wave table
    background <- original_rate_df %>%
      dplyr::filter(group=="Background") %>%
      dplyr::group_by(measurement) %>%
      dplyr::reframe(bkg_OCR_wave = mean(OCR_wave),
                       bkg_ECAR_wave = mean(ECAR_wave)
      )
    original_rate_df <- dplyr::left_join(original_rate_df,
                                         background,
                                         by = c("measurement"), copy = TRUE)

    original_rate_df$OCR_wave_bc <- original_rate_df$OCR_wave - original_rate_df$bkg_OCR_wave
    original_rate_df$ECAR_wave_bc <- original_rate_df$ECAR_wave - original_rate_df$bkg_ECAR_wave

    original_rate_df <- original_rate_df %>%
      dplyr::select(measurement, well, group,
                    time_wave, OCR_wave, OCR_wave_bc,
                    ECAR_wave, ECAR_wave_bc)
  }

  original_rate_df_list <- list(original_rate_df, corrected_allready)

  logger::log_info("Finished collecting OCR from 'Rate' sheet.")

  return(original_rate_df_list)

}

check_excel_positions <- function(df, pos_vector, name_vector){

  logger::log_info("Check if excel df contains data name on certain position.")
  tf_values <- mapply(function(pos, name) {
    true_false <- name %in% df[[1]][pos]
    return(true_false)
  }, pos_vector, name_vector)

  check_tf_list <- function(tf_values){
    if (all(tf_values)) {
      logger::log_info("'Assay Configuration' sheet contains all values.")
      return(TRUE)
    } else {
      logger::log_error("'Assay Configuration' sheet doesn't contain all values.")
      stop()
    }
  }

  tf <- check_tf_list(tf_values)

  return(tf)
}

get_xf_rate <- function(filepath_seahorse){
    #first item is table, second item is background_corrected logical
    xf_rate_list <- get_originalRateTable(filepath_seahorse)

    if (exists("xf_rate_list")) {
      cli::cli_alert_success("Succesfully read Rate data from data sheet.")
    } else {
      stop("An error occured during the data collection from 'Rate'sheet.")
    }

    logger::log_info("Finished collecting data from 'Rate' sheet.")

    return(xf_rate_list)
}

get_xf_buffer <- function(filepath_seahorse){

    logger::log_info("Collecting buffer factor info from 'Assay configuration' sheet.")


    bufferfactor_info <- get_platelayout_data(filepath_seahorse,
                                              my_sheet = "Assay Configuration",
                                              my_range = "B96:N104",
                                              my_param = "bufferfactor")

    logger::log_info("Finished collecting buffer factor info from 'Assay configuration' sheet.")

    if (exists("bufferfactor_info")) {
      cli::cli_alert_success("Succesfully read bufferfactor information from 'Assay configuration' data sheet.")
    } else {
      stop("An error occured during the data collection from 'Assay Configuration' sheet.")
    }

    logger::log_info("Finished collecting bufferfactor information from 'Assay configuration' sheet.")


    return(bufferfactor_info)

}

get_xf_pHcal <- function(filepath_seahorse){
  logger::log_info("Collecting pH calibration emission data")

  pH_calibration <- get_platelayout_data(filepath_seahorse,
                                         my_sheet = "Calibration",
                                         my_range = "P16:AB24",
                                         my_param = "pH_cal_em")



  if (exists("pH_calibration")) {
    cli::cli_alert_success("Succesfully read pH calibration emission data from 'Calibration' data sheet.")
  } else {
    stop("An error occured during the data collection from 'Calibration sheet.")
  }

  logger::log_info("Finished collecting pH calibration emission data")


  return(pH_calibration)
}

get_xf_O2cal <- function(filepath_seahorse){

  logger::log_info("Collecting O2 calibration emission")

  O2_calibration <- get_platelayout_data(filepath_seahorse,
                                         my_sheet = "Calibration",
                                         my_range = "B7:N15",
                                         my_param = "O2_cal_em")

  if (exists("O2_calibration")) {
    cli::cli_alert_success("Succesfully collected O2 calibration data.")
  } else {
    stop("An error occured during the data collection from 'Calibration sheet.")
  }

  logger::log_info("Finished collecting O2 calibration emission")


  return(O2_calibration)
}

get_xf_inj <- function(filepath_seahorse, injscheme = "HAP"){

  logger::log_info("Collecting injection information")

  #command_index in "Operation Log" sheet give numbers to the phases in a seahorse exp
  # each command (eg. "mix", "measure") gets the command_index for that phase
  # 0 = moving operation
  # 1 = calibration
  # 2 = baseline
  # 3 = injection 1
  # 4 = injection 2
  # 5 = injection 3
  # 6 = injection 4

  #read injection strategy and measurements from "Operation Log" sheet
  info_sh<-readxl::read_excel(filepath_seahorse, sheet = "Operation Log")
  colnames(info_sh) <- c("instruction_name","command_name",
                         "command_index","start_time","end_time",
                         "completion_status")

  if (injscheme == "HAP"){
    #assumes injection names are available in operation log file (this is the case for most experiments)
    measurement_info <- dplyr::filter(info_sh, command_name == "Measure")
    measurement_info$interval <- measurement_info$command_index -1
    measurement_info$measurement <- 1:nrow(measurement_info)
    measurement_info <- measurement_info %>% dplyr::select(measurement, interval, injection=instruction_name)
  }

  if (injscheme == "manual"){

    #in case there is no command index in "opertion log"
    command_names <- c("XF - PC_Measure", "XF - PC_Inject")
    measurement_info <- dplyr::filter(info_sh, command_name %in% command_names)

    # "PC - inject" has a number as command_index
    # "PC - measure" command_index == 0
    # I use that to set the command_index
    interval = 1
    for (i in 1:nrow(measurement_info)){
      if(measurement_info$command_index[i] == 0){
        measurement_info$command_index[i] <-  interval } else {
          interval <-  interval +1
          measurement_info$command_index[i] <-  interval}
    }
    colnames(measurement_info)[3] <- "interval"
    measurement_info <- dplyr::filter(measurement_info, command_name == "XF - PC_Measure")
    measurement_info$measurement <- 1:nrow(measurement_info)
    measurement_info <- measurement_info %>% dplyr::select(measurement, interval)

    #gives name of the injection manually
    # case mitostress
    injections <- c("basal", "OM", "FCCP", "AM/rot")
    injections_mitostress <- tibble::tibble(interval = 1:4, injection=c("basal", "OM", "FCCP", "AM/rot"))
    measurement_info <- dplyr::left_join(measurement_info, injections_mitostress, by = c("interval"))

    ## case glycostress
    #injections <- c("basal", "glucose", "OM", "2DG")
    #injections_glycostress <- tibble(interval = 1:4, injection=injections)
    #measurement_info <- left_join(measurement_info, injections_glycostress, by = c("interval"))
  }

  logger::log_info("Finished collecting injection information")


  if (exists("measurement_info")) {
    cli::cli_alert_success("Succesfully read injection information from 'Assay configuration' data sheet.")
  } else {
    stop("An error occured during the data collection from 'Assay Configuration' sheet.")
  }

  logger::log_info("Finished collecting injection information  from 'Assay configuration' sheet.")

  return(measurement_info)

}

get_xf_assayinfo <- function(filepath_seahorse,
                             date_style = "empty",
                             instrument = "XFe96",
                             norm_available,
                             xls_ocr_backgroundcorrected) {

  logger::log_info("Collecting assay information")

  if (instrument == "XFHSmini"){
    gain1_cell <-  "D68"
    gain2_cell <-  "E68"
  }

  if (instrument == "XFe96"){
    gain1_cell <-  "D70"
    gain2_cell <-  "E70"
  }

  # read Assay Configuration sheet
  try(
    meta_df <- readxl::read_excel(filepath_seahorse,
                          sheet = "Assay Configuration",
                          col_names = c("parameter", "value"),
                          range = "A1:B83"
    )
  )

  if (!(exists("meta_df"))) {
    stop("An error occured during the data collection from 'Assay Configuration' sheet.")
  }

  pos_vector = c(4, 26, 32, 38, 58, 59, 60, 61, 62, 63, 65, 66, 67, 76)
  name_vector = c("Assay Name", "Cartridge Barcode", "Plate Barcode", "Instrument Serial",
                  "ksv", "Ksv Temp Correction", "Corrected Ksv", "Calculated FO",
                  "Pseudo Volume", "TAC", "TW", "TC", "TP", "Calibration pH")

  # Assertion to check if meta_df contains Seahorse constants on the right locations.
  tf <- check_excel_positions(meta_df, pos_vector, name_vector)

  meta_df <- meta_df[!is.na(meta_df$parameter), ]

  # read Assay Configuration sheet gain1
  try(
    gain1 <- readxl::read_excel(filepath_seahorse,
                        sheet = "Assay Configuration",
                        col_names = c("value"),
                        range = gain1_cell
    )
  )

  if (!(exists("gain1"))) {
    stop("An error occured during the data collection from 'Assay Configuration' sheet.")
  }


  # read Assay Configuration sheet gain2
  try(
    gain2 <- readxl::read_excel(filepath_seahorse,
                        sheet = "Assay Configuration",
                        col_names = c("value"),
                        range = gain2_cell
    )
  )
  if (!(exists("gain2"))) {
    stop("An error occured during the data collection from 'Assay Configuration' sheet.")
  }



  # read target emission cells
  try(
  O2_target_emission <- readxl::read_excel(filepath_seahorse,
                                   sheet = "Calibration",
                                   col_names = FALSE,
                                   range = "B4"
  ))

  if (!(exists("O2_target_emission"))) {
    stop("An error occured during the data collection from 'Calibration' sheet.")
  }


  # read pH target emission cells
  try(
    pH_target_emission <- readxl::read_excel(filepath_seahorse,
                                     sheet = "Calibration",
                                     col_names = FALSE,
                                     range = "P4"


    ))
  if (!(exists("gain1"))) {
    stop("An error occured during the data collection from 'Calibration' sheet.")
  }

  F0 <- as.double(meta_df$value[meta_df$parameter == "Calculated FO"])
  V_C <- as.double(meta_df$value[meta_df$parameter == "Pseudo Volume"])
  Tau_AC <- as.double(meta_df$value[meta_df$parameter == "TAC"])
  Tau_W <- as.double(meta_df$value[meta_df$parameter == "TW"])
  Tau_C <- as.double(meta_df$value[meta_df$parameter == "TC"])
  Tau_P <- as.double(meta_df$value[meta_df$parameter == "TP"])
  KSV_original <- as.double(meta_df$value[meta_df$parameter == "ksv"])
  KSV_corrected <- as.double(meta_df$value[meta_df$parameter == "Corrected Ksv"])
  KSV_tempCorrection <- as.logical(meta_df$value[meta_df$parameter == "Ksv Temp Correction"])
  KSV <- KSV_corrected

  pH_0 <- as.double(meta_df$value[meta_df$parameter == "Calibration pH"])
  pH_plateVolume <- as.double(meta_df$value[meta_df$parameter == "Plate Volume"])
  pH_kVol <- as.double(meta_df$value[meta_df$parameter == "kVol"])


  plate_id <- meta_df$value[meta_df$parameter == "Plate Barcode"]
  cartridge_barcode <- meta_df$value[meta_df$parameter == "Cartridge Barcode"]
  assay_name <- meta_df$value[meta_df$parameter == "Assay Name"]
  instrument_serial <- meta_df$value[meta_df$parameter == "Instrument Serial"]

  pH_targetEmission <- as.double(pH_target_emission[[1]])
  O2_targetEmission <- as.double(O2_target_emission[[1]])
  gain1 <- as.double(gain1[[1]])
  gain2 <- as.double(gain2[[1]])

  # other constants
  O2_0_mmHg <- 151.6900241
  O2_0_mM <- 0.214

  if (date_style == "US"){
    date_run <- lubridate::mdy_hm(meta_df$value[meta_df$parameter == "Last Run"])
    logger::log_info("Converted date to US format (US = mdy_hm, NL = dmy_hm).") # (Date-time column)
    #be carefull with the data format in excel! either mdy or dmy
  }

  if (date_style == "NL"){
    date_run <- lubridate::dmy_hm(meta_df$value[meta_df$parameter == "Last Run"])
    logger::log_info("Converted date to NL format (US = mdy_hm, NL = dmy_hm).") # (Date-time column)
    #be carefull with the data format in excel! either mdy or dmy
  }

  if (date_style == "empty"){
    date_run <- meta_df$value[meta_df$parameter == "Last Run"] # (Character instead of date-time column)
    logger::log_info("Date-style is empty, no date conversion was performed. Format is 'character' instead of 'date'.")
    #be carefull with the data format in excel! either mdy or dmy
  }

  if(instrument == "XFHSmini"){
    tibbler <- tibble::tibble(
      F0 = 4.63e04,
      V_C = 9.15,
      Tau_AC = 746,
      Tau_W = 296,
      Tau_C = 246,
      Tau_P = 60.9,
      KSV = 2.06e-02,
      KSV_corrected = 2.06e-02,
      KSV_original = 2.06e-02,
      KSV_tempCorrection = FALSE,
      gain1,
      gain2,
      pH_0,
      pH_plateVolume,
      pH_kVol,
      pH_targetEmission,
      O2_targetEmission,
      plate_id,
      cartridge_barcode,
      date_run,
      assay_name,
      instrument_serial,
      O2_0_mmHg,
      O2_0_mM
    )

  }
  if(instrument == "XFe96"){
    tibbler <- tibble::tibble(
      F0,
      V_C,
      Tau_AC, Tau_W,
      Tau_C, Tau_P,
      KSV,
      KSV_tempCorrection,
      KSV_original,
      gain1,
      gain2,
      pH_0,
      pH_plateVolume,
      pH_kVol,
      pH_targetEmission,
      O2_targetEmission,
      plate_id,
      cartridge_barcode,
      date_run,
      assay_name,
      instrument_serial,
      O2_0_mmHg,
      O2_0_mM
    )
  }


  tibbler$norm_available <- norm_available
  tibbler$xls_ocr_backgroundcorrected <- xls_ocr_backgroundcorrected

  cli::cli_progress_bar("Finished collecting assay information.")

  return(tibbler)
}

get_platelayout_data <- function(filepath_seahorse, my_sheet, my_range, my_param ){

      df <- readxl::read_excel(filepath_seahorse, sheet = my_sheet, range = my_range)

      colnames(df)[1] <- "firstCol"

      df <-  tidyr::gather(df, key = "key", value = "my_value", -firstCol) %>%
        dplyr::mutate(firstCol = paste0(firstCol, key) ) %>%
        dplyr::select(well = firstCol, my_value) %>%
        dplyr::arrange(gsub("\\d", "", well, as.numeric(gsub("\\D", "", well))))

      colnames(df)[2] <- my_param

      # add a zero between letter and number if wellname has 2 characters for normalization data
      for (i in 1:nrow(df)){
        if (nchar(df$well[i]) ==  2) {
          wellName <- sub("(.{1})(.*)", "\\10\\2", df$well[i])
        } else {
          wellName <- df$well[i]
        }
        df$well[i] <- wellName
      }

     return(df)

}



