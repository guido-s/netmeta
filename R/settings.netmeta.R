#' Print and change default network meta-analysis settings in R package
#' \bold{netmeta}
#' 
#' @description
#' Print and change default settings to conduct and print or plot
#' network meta-analyses in R package \bold{netmeta}.
#' 
#' @param ... Arguments to change default settings.
#' @param quietly A logical indicating whether information on settings
#'   should be printed.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{settings.meta}}, \code{\link[meta]{gs}}
#' 
#' @export settings.netmeta

settings.netmeta <- function(..., quietly = TRUE) {
  
  #
  # Check argument
  #
  
  missing.quietly <- missing(quietly)
  chklogical(quietly)
  
  
  #
  # Save object with current settings
  #
  
  argslist <- gs(".argslist.netmeta")
  #
  oldset <- vector("list", 0)
  #
  for (i in argslist)
    oldset[[i]] <- gs(i)
  #
  print.settings <- FALSE
  reset.settings <- FALSE
  #
  args  <- list(...)
  
  names.all <- names(args)
  #
  sel.netmeta <- names.all %in% argslist
  
  
  #
  # Determine whether to print, reset or use specific settings
  #
  if (any(names.all == "print") && args[["print"]]) {
    print.settings <- TRUE
    quietly <- FALSE
  }
  if (any(names.all == "reset") && args[["reset"]])
    reset.settings <- TRUE
  
  
  #
  # Print settings if no argument is provided
  #
  
  if (length(args) == 0) {
    if (missing.quietly || !quietly)
      settings.netmeta("print", quietly = FALSE)
    #
    return(invisible(oldset))
  }
  
  
  #
  # Check argument names
  #
  
  names.all <- names(args)
  #
  if (length(names.all) != length(unique(names.all)))
    stop("Arguments must be unique.")
  #
  sel.netmeta <- names.all %in% argslist
  names.netmeta <- names.all[sel.netmeta]
  names.other <- names.all[!sel.netmeta]
  #
  args.netmeta <- args.other <- vector("list", 0)
  #
  if (length(names.netmeta) > 0) {
    for (i in names.netmeta) {
      args.netmeta[[i]] <- args[[i]]
    }
  }
  #
  if (length(names.other) > 0) {
    for (i in names.other) {
      args.other[[i]] <- args[[i]]
    }
  }
  
  
  #
  # Check whether first argument is a list. In this case only use
  # this list as input.
  #
  
  warn.depr <- TRUE
  if (length(args) > 0 && is.list(args[[1]])) {
    if (!is.null(names(args))) {
      print(names(args))
      warning("Additional arguments ignored as first argument is a list.",
              call. = FALSE)
    }
    warn.depr <- FALSE
    args <- args[[1]]
  }
  
  
  #
  # Unnamed first (and only) argument must be character string or a
  # logical
  #
  
  if (length(args) == 1 & is.null(names(args))) {
    if (is.character(unlist(args)))
      action <- setchar(unlist(args), c("reset", "print"),
                        stop.at.error = FALSE)
    else
      action <- unlist(args)
    #
    if (is.null(action))
      stop("First argument can be one of the following character strings:",
           "\n 'reset', 'print'",
           call. = FALSE)
    else if (action == "reset")
      settings.netmeta(reset = TRUE, quietly = quietly)
    else if (action == "print" | (is.logical(action) && action))
      settings.netmeta(print = TRUE, quietly = FALSE)
    #
    return(invisible(oldset))
  }
  #
  else if (length(args) > 1 & names(args)[1] == "") {
    if (is.character(unlist(args[[1]])))
      action <- setchar(unlist(args[[1]]), c("reset", "print"),
                        stop.at.error = FALSE)
    else
      action <- unlist(args[[1]])
    #
    if (is.null(action))
      stop("First argument can be one of the following character strings:",
           "\n 'reset', 'print'",
           call. = FALSE)
    else if (action == "reset")
      settings.netmeta(reset = TRUE, quietly = quietly)
    else if (action == "print")
      settings.netmeta(print = TRUE, quietly = FALSE)
  }
    
  
  #
  # Reset settings
  #
  
  if (reset.settings) {
    if (!quietly)
      cat("\n** Reset all network meta-analysis settings (R package netmeta). ",
          "**\n\n")
    #
    settings.meta(baseline.reference = TRUE, small.values = "desirable",
                  all.treatments = NULL, seq = NULL,
                  method.tau.netmeta = "DL",
                  drop.reference.group = TRUE, equal.size = TRUE,
                  show = "both",
                  #
                  nsim = 1000, lump.comparator = FALSE,
                  #
                  plastic = FALSE, col.netgraph = NULL,
                  number.of.studies = TRUE, thickness = "number.of.studies",
                  multiarm = FALSE,
                  #
                  tol.multiarm = 0.001, tol.multiarm.se = NULL,
                  details.chkmultiarm = FALSE,
                  #
                  na.unident = TRUE,
                  sep.trts = ":", sep.comps = "+", sep.ia = "*",
                  nchar.trts = 666, nchar.studlab = 666,
                  #
                  legend = TRUE)
  }
  
  
  #
  # Print settings
  #
  
  if (print.settings & !quietly) {
    cat(paste0("\n** Settings for network meta-analysis method (",
               "R package netmeta, version ",
               utils::packageDescription("netmeta")$Version,
               ") **\n\n"))
    #
    cat(paste0("* General settings *\n"))
    catarg("baseline.reference   ")
    catarg("small.values         ")
    catarg("all.treatments       ")
    catarg("seq                  ")
    catarg("method.tau.netmeta   ")
    catarg("tol.multiarm         ")
    catarg("tol.multiarm.se      ")
    catarg("details.chkmultiarm  ")
    catarg("nchar.studlab        ")
    catarg("legend               ")
    #
    cat("\n* Additional settings for network meta-analysis *\n")
    catarg("sep.trts             ")
    catarg("nchar.trts           ")
    #
    cat("\n* Additional settings for component network meta-analysis *\n")
    catarg("sep.comps            ")
    catarg("sep.ia               ")
    catarg("na.unident           ")
    #
    cat("\n* Additional settings for forest plots *\n")
    catarg("drop.reference.group ")
    catarg("equal.size           ")
    #
    cat("\n* Additional settings for network graphs *\n")
    catarg("plastic              ")
    catarg("col.netgraph         ")
    catarg("number.of.studies    ")
    catarg("thickness            ")
    catarg("multiarm             ")
    #
    cat("\n* Additional setting for netsplit() *\n")
    catarg("show                 ")
    #
    cat("\n* Additional setting to evaluate small study effects *\n")
    catarg("lump.comparator      ")
    #
    cat("\n* Number of samples to calculate ranking metrics *\n")
    catarg("nsim                 ")
    #
    return(invisible(oldset))
  }
  
  
  #
  # Set settings
  #

  if (length(names.netmeta) > 0) {
    setlogical("baseline.reference", args.netmeta)
    setcharacter("small.values", args.netmeta, c("desirable", "undesirable"))
    setlogical("all.treatments", args.netmeta, NULL.ok = TRUE)
    setcharacter("seq", args.netmeta, NULL.ok = TRUE)
    setcharacter("method.tau.netmeta", args.netmeta, c("DL", "REML", "ML"))
    setlogical("drop.reference.group", args.netmeta)
    setlogical("equal.size", args.netmeta)
    setcharacter("show", args.netmeta,
                 c("all", "both", "with.direct", "direct.only", "indirect.only"))
    setnumeric("nsim", args.netmeta)
    setlogical("lump.comparator", args.netmeta)
    setlogical("plastic", args.netmeta)
    setcharacter("col.netgraph", args.netmeta, NULL.ok = TRUE)
    setlogical("number.of.studies", args.netmeta)
    setcharacter("thickness", args.netmeta)
    setlogical("multiarm", args.netmeta)
    setnumeric("tol.multiarm", args.netmeta)
    setnumeric("tol.multiarm.se", args.netmeta)
    setlogical("details.chkmultiarm", args.netmeta)
    setlogical("na.unident", args.netmeta)
    setcharacter("sep.trts", args.netmeta)
    setcharacter("sep.comps", args.netmeta)
    setcharacter("sep.ia", args.netmeta)
    setnumeric("nchar.trts", args.netmeta)
    setnumeric("nchar.studlab", args.netmeta)
    setlogical("legend", args.netmeta)
  }
  
  
  #
  # Return current settings
  #
  
  res <- vector("list", 0)
  #
  for (i in argslist)
    res[[i]] <- gs(i)
  #
  invisible(res)
}
