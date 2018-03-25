
str_to_lang <- function(str) {
    if (is.character(str)) {
        return(parse(text = str)[[1]])
    } else {
        return(str)
    }
}

symbol_list <- function(lang) {
    res <- lapply(lang, function(el) {
        if (!is.language(el)) {
            # drop constants
            NULL
        } else if (is.symbol(el)) {
            # convert symbol to string
            as.character(el)
        } else {
            # have language that's still more than one symbol, so recurse
            symbol_list(el)
        }
    })
    unlist(res)
}


check_symbol_list <- function(symbols, data = list()) {
    data_env <- list2env(data, parent = parent.env(environment()))
    check <- sapply(symbols, exists, envir = data_env)
    missing <- symbols[!check]
    if (length(missing) > 0) {
        warning(paste("Missing symbols: ", paste(missing, collapse = ", ")))
        return(FALSE)
    } else {
        return(TRUE)
    }
}

# voodoo to partially evaluate an expression
expr_sub <- function(expr, env) {
    subbed <- sapply(expr, function(e1) {
        eval(substitute(substitute(e, env), list(e = e1)))
    })
    as.expression(subbed)
}

define_param <- function(psi, HA, CA, IC) {
    # get default
    args <- formals()

    # replace defaults with user arguments when specified
    user_args <- as.list(match.call(expand.dots = TRUE))[-1]
    args[names(user_args)] <- user_args

    # verify we have language objects, not strings
    sapply(args, str_to_lang)
}
