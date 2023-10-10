# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {

}

n <- 1000
data <- tibble::tibble(
  X2 = as.factor(sample(1:2,size = n,TRUE)),
  X4 = as.character(sample(1:4,size = n,TRUE)),
  aa = rnorm(n, mean=10,sd=1),
  bb = rnorm(n, mean=10,sd=1),
  cc = rnorm(n, mean=10,sd=1),
  dd = rnorm(n, mean=10,sd=1),
  ee = as.character(sample(1:2,size = n,TRUE)),
  ff = as.factor(sample(1:4,size = n,TRUE)),
  gg = sample(1:3,size = n,TRUE)
)

#data <- tibble::as_tibble(datasets::iris)

#--args--------------------
roworder <- c("aa","bb","dd","ee", "ff", "gg")
byvars <- c("X2")
stats <- list("{mean}±{sd}" = c("aa","dd"),
              "{p50}[{p25}-{p75}]" = c("bb")
              )
stats_accuracy <- list(default = 0.01, aa = 0.001)

counts <- list("{n}({perc})" = c("ee","ff"),
               "{n}/{total}" = c("gg"))
percent_accuracy <- 0.1

output = "md" #or "flextable"
#--------------------------

#利用する関数群（隠す）
make_single_column <- function(data, stats, stats_accuracy, counts, percent_accuracy){
  numeric_data <- purrr::map_dfr(1:length(stats), ~{
    gluestring <-names(stats[.])
    targetvars <-stats[[.]]

    funs_in_gluestring <- stringr::str_extract_all(gluestring, "\\{mean\\}|\\{sd\\}|\\{p\\d+\\}|\\{min\\}|\\{max\\}")[[1]] |>
      stringr::str_remove_all("\\{|\\}")

    individual_data <- purrr::map(funs_in_gluestring, ~{
      afunstr <- .

      if(afunstr == "mean") fun <- mean
      if(afunstr == "sd") fun <- sd
      if(afunstr == "min") fun <- min
      if(afunstr == "max") fun <- max
      if(stringr::str_detect(afunstr,"p\\d+")) {
        num <- as.numeric(stringr::str_extract(afunstr,"\\d+"))
        fun <- function(x){
          quantile(x, p = num/100)
        }
      }

      fin <- data |>
        dplyr::summarise(dplyr::across(all_of(targetvars), fun, .names = "{.col}"))

      fin <- fin |>
        tidyr::pivot_longer(cols = everything(), names_to = "col") |>
        dplyr::mutate(value = purrr::map2_chr(col, value, ~{
          accu <- stats_accuracy$default
          if(.x %in% names(stats_accuracy)) accu <- stats_accuracy[[.x]]
          return(scales::number(.y, accu))
        })) |>
        dplyr::rename(!!rlang::sym(afunstr) := value)


      return(fin)

    })

    if(length(individual_data) > 1){
      combined_data <- individual_data |> purrr::reduce(~{
        dplyr::left_join(.x,.y,by="col")
      })
    }else{
      combined_data <- individual_data
    }

    combined_data <- combined_data |>
      dplyr::mutate(res = stringr::str_glue(gluestring))

    return(combined_data)
  })
  numeric_data <- numeric_data |>
    dplyr::mutate(type = NA_character_) |>
    dplyr::select(var = col, type, res)


  count_data <- purrr::map_dfr(1:length(counts), ~{
    gluestring <-names(counts[.])
    targetvars <-counts[[.]]

    funs_in_gluestring <- stringr::str_extract_all(gluestring, "\\{n\\}|\\{perc\\}|\\{total\\}")[[1]] |>
      stringr::str_remove_all("\\{|\\}")

    individual_data <- purrr::map_dfr(targetvars, ~{

      atgt <- .
      data |>
        dplyr::count(!!rlang::sym(atgt)) |>
        dplyr::add_row(.before = 1) |>
        dplyr::mutate(var = atgt) |>
        dplyr::rename(type = !!rlang::sym(atgt)) |>
        dplyr::mutate(type = as.character(type)) |>
        dplyr::relocate(var,type,n) |>
        dplyr::mutate(total = sum(n, na.rm=TRUE),
                      perc = scales::number(100*n/total,accuracy=percent_accuracy)) |>
        dplyr::mutate(res = stringr::str_glue(gluestring)) |>
        dplyr::select(var, type, res) |>
        dplyr::mutate(res = dplyr::if_else(is.na(type),NA_character_, as.character(res)))

    })

    return(individual_data)
  })

  single_column_data <- dplyr::bind_rows(numeric_data, count_data)

  return(single_column_data)
}
make_base_table <- function(data,  byvars, stats, stats_accuracy, counts, percent_accuracy){
  browser()
  total_column <- make_single_column(data, stats, stats_accuracy, counts, percent_accuracy)

  totaldata <- total_column |>
    dplyr::mutate(header = "total") |>
    dplyr::group_by(header) |>
    tidyr::nest() |>
    dplyr::select(header, columns = data) |>
    dplyr::mutate(n = nrow(data))

  if(length(byvars) == 1){
    if(byvars == ""){
      # do nothing
      stopifnot("column name , '...temporary_group_zero2346598' is used in the package please rename in prior to use epitablex package" = !("...temporary_group_zero2346598" %in% colnames(data)))
      grpdata <- data |> dplyr::mutate(...temporary_group_zero2346598 = "...temporary_group_zero2346598") |> dplyr::group_nest(...temporary_group_zero2346598)
      byvars <- "...temporary_group_zero2346598"
    }else{
      grpdata <- data |> dplyr::group_nest(!!!rlang::syms(byvars))
    }
  }else{
    grpdata <- data |> dplyr::group_nest(!!!rlang::syms(byvars))
  }

  grpdata <- grpdata |>
    dplyr::mutate(columns = purrr::map(data, ~{
      adata <- .
      make_single_column(adata , stats, stats_accuracy, counts, percent_accuracy)
    })) |>
    dplyr::mutate(n = purrr::map_int(data,nrow)) |>
    dplyr::select(!data)

  grpdata <- grpdata |>
    dplyr::rowwise() |>
    dplyr::mutate(header = stringr::str_c(dplyr::c_across(byvars), collapse="_")) |>
    dplyr::select(header, columns, n)

  grpdata <- grpdata |> dplyr::ungroup()

  finaldata <- dplyr::bind_rows(totaldata, grpdata)

  if("...temporary_group_zero2346598" %in% finaldata$header){
    finaldata <- finaldata |> dplyr::filter(header != "...temporary_group_zero2346598")
  }

  return(finaldata)
}
make_wide_table <- function(basetable, add_n = "n = {n}", pos_total = "left", roworder = NA){

  stopifnot("pos_total need to be left, right or none" = pos_total %in% c("left","right","none"))

  #make widetable from basetable data.
  widetable <- purrr::map2(basetable$header, basetable$columns, ~{
    dplyr::rename(.y, !!rlang::sym(.x) := res)
  }) |>
    purrr::reduce(~{
      dplyr::left_join(.x, .y, by=c("var","type"))
    })

  #change order of row variables
  if(any(is.na(roworder))){
    #do nothing
  }else{
    stopifnot("All variable present in arguments stats and counts are need to be present in argument roworder." = all(unique(widetable$var) %in% roworder) )

    widetable <- widetable |>
      dplyr::mutate(var = factor(var, levels = roworder)) |>
      dplyr::arrange(var) |>
      dplyr::mutate(var = as.character(var))
  }

  #Change position of total or omit total
  if(pos_total == "left"){
    widetable <- widetable |>  dplyr::relocate(var, type, total)
  }else if(pos_total == "right"){
    widetable <- widetable |> dplyr::relocate(total, .after = dplyr::last_col())
  }else if(pos_total == "none"){
    widetable <- widetable |> dplyr::select(!total)
  }

  #tidy widetable data
  widetable <- widetable |>
    dplyr::mutate(var = dplyr::if_else(is.na(type),var, as.character(stringr::str_glue("\t{type}")))) |>
    dplyr::select(!type)


  return(widetable)
}

gen_nrow <- function(basetable, widetable, add_n = "N = {n}"){
  nvalues <- basetable |>
    dplyr::select(header, n) |>
    dplyr::filter(header %in% colnames(widetable)) |>
    dplyr::mutate(val = stringr::str_glue(add_n)) |>
    dplyr::select(header, val) |>
    dplyr::ungroup() |>
    dplyr::add_row(tibble::tibble(header = "var", val = "")) |>
    dplyr::mutate(header = factor(header, levels = colnames(widetable))) |>
    dplyr::arrange(header)

  return(nvalues)
}

make_output_flextable <- function(basetable, widetable, add_n = "N = {n}"){

  ft <- flextable::flextable(widetable) |>
    flextable::autofit()

  if(!is.character(add_n)){
    #if add_n is FALSE then do nothing
  }else{
    nvalues <- gen_nrow(basetable, widetable, add_n)
    ft <- ft |>flextable::add_header_row(nvalues$val, top=FALSE)
  }

  return(ft)

}

make_output <- function(basetable, widetable, add_n = "N = {n}", as_md = FALSE){
  browser()
  result <- widetable

  if(!is.character(add_n)){
    #if add_n is FALSE then do nothing
  }else{
    nvalues <- gen_nrow(basetable, widetable, add_n)
    rowtibble <- nvalues |>
      dplyr::mutate(..id = 1) |>
      tidyr::pivot_wider(id_cols = ..id, names_from = header, values_from = val) |>
      dplyr::select(!..id)

    result <- result |> dplyr::add_row(rowtibble, .before=1)
  }

  if(as_md){
    result <- result |> knitr::kable()
  }

  return(result)
}

epitable1 <- function(data,
                      byvars= "",
                      roworder = "",
                      stats = "{mean}±{sd}",
                      stats_accuracy = 0.01,
                      counts = "{n}({perc})",
                      percent_accuracy = 0.1,
                      add_n = "n = {n}", pos_total = "left", output = "md"){

  default_stat  <- "{mean}±{sd}"
  default_count <- "{n}({perc})"


  #PROCESS ROWORDER--------------------------------------------------
  #if roworder=="" then make roworder depend on colnames of data
  if(all(roworder %in% "")){
    roworder <- colnames(data)
    roworder <- roworder[!roworder %in% byvars]

    if(length(byvars) > 0) roworder <- roworder[!roworder %in% byvars]
  }else{
    if(all(byvars %in% "")){
      selection <- c(roworder)
    }else{
      selection <- c(roworder, byvars)
    }

    data <- data |> dplyr::select(!!!rlang::syms(selection))
  }

  #PROCESS STATS-------------------------------------------------------
  #if non list is supplied to stats, then make list type stats arguments
  #stats is applied to numeric or integer columns.
  if(class(stats)=="character"){
    tgts <- data |> dplyr::select(where(~{is.numeric(.) | is.integer(.) })) |> colnames()
    if(length(byvars) > 0) tgts <- tgts[!tgts %in% byvars]
    stats <- setNames(list(c(tgts)), stats)
  }else if(class(stats)=="list"){
    #if stats and roworder have different number of variables specified, then set other to
    #default stat "mean(sd)" value
    tgts <- data |> dplyr::select(where(~{is.numeric(.) | is.integer(.) })) |> colnames()
    notspecified <- tgts[!tgts %in% purrr::flatten_chr(stats)]
    if(length(byvars) > 0) notspecified <- notspecified[!notspecified %in% byvars]
    if(default_stat %in% names(stats)){
      stats <- append(stats, setNames(list(notspecified), default_stat))
    }else{
      stats[[default_stat]] <- c(stats[[default_stat]], notspecified)
    }
  }

  #PROCESS STATS ACCURACY----------------------------------------------
  #if only single value supplied, then set that as default
  if(is.numeric(stats_accuracy) & length(stats_accuracy)==1){
    stats_accuracy <- list(default = stats_accuracy)
  }

  #PROCESS COUNTS-------------------------------------------------------
  #if non list is supplied to counts, then make list type counts arguments
  #counts is applied to character or factor columns.
  if(class(counts)=="character"){
    tgts <- data |> dplyr::select(where(~{is.factor(.) | is.character(.) })) |> colnames()
    if(length(byvars) > 0) tgts <- tgts[!tgts %in% byvars]
    counts <- setNames(list(c(tgts)), counts)
  }else if(class(counts)=="list"){
    #if counts and roworder have different variables specified, then set other to
    tgts <- data |> dplyr::select(where(~{is.factor(.) | is.integer(.) })) |> colnames()
    notspecified <- tgts[!tgts %in% purrr::flatten_chr(counts)]
    if(length(byvars) > 0) notspecified <- notspecified[!notspecified %in% byvars]
    if("{n}({perc})" %in% names(counts)){
      counts <- append(counts, setNames(list(notspecified), default_count))
    }else{
      counts[[default_count]] <- c(counts[[default_count]], notspecified)
    }
  }

  #PROCESS byvars-----------------------------------------------------------
  if(length(byvars) == 1){
    if(byvars == ""){
      #do nothing
    }else{
      stopifnot("byvars should be colnames in data" = byvars %in% colnames(data))
    }
  }else{
    stopifnot("byvars should be colnames in data" = all(byvars) %in% colnames(data))
  }

  basetable <- make_base_table(data,  byvars, stats, stats_accuracy, counts, percent_accuracy)
  widetable <- make_wide_table(basetable, add_n = add_n, pos_total = pos_total, roworder = roworder)

  if(output == "flextable") result <- make_output_flextable(basetable, widetable, add_n)
  if(output == "md") result <- make_output(basetable, widetable, add_n, as_md=TRUE)
  if(output == "tibble") result <- make_output(basetable, widetable, add_n, as_md=FALSE)

  return(result)
}


data |>
  #dplyr::select(!X4) |>
  epitable1(
    byvars = "X4",
    stats = list("{mean}±{sd}" = c("aa","bb"),
                 "{p50}[{p25}-{p75}]" = c("cc","dd")),
    pos_total = "right",
    output="md")





roworder <- c("aa","bb","dd","ee", "ff", "gg")
byvars <- c("X2")
stats <- list("{mean}±{sd}" = c("aa","dd"),
              "{p50}[{p25}-{p75}]" = c("bb")
)
stats_accuracy <- list(default = 0.01, aa = 0.001)

counts <- list("{n}({perc})" = c("ee","ff"),
               "{n}/{total}" = c("gg"))
percent_accuracy <- 0.1

output = "md" #or "flextable"


