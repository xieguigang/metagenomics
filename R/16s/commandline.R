#' the mothur cli wrapper
#' 
#' @details Schloss PD et al. 2009. Introducing mothur: Open-source, 
#'   platform-independent, community-supported software for describing 
#'   and comparing microbial communities. Applied and Environmental 
#'   Microbiology 75:7537–7541. 
#'
#' @param command the mothur command string.
#' @param log the logfile file path.
#'
const runMothur as function(command, argv, log) {
    const mothur as string      = getOption("mothur");
    const is_debug as boolean   = getOption("workflow.debug");
    const commandArgs as string = `${mothur} "#${command}(${command |> mothurArgvs(argv)})"`;

    cat(`${log} ${commandArgs}\n`);

    if (!is_debug) {
        commandArgs
        |> system(intern = TRUE)
        |> writeLines(con = log)
        ;

        stdout = readLines(log);

        if (any(stdout like $".*\[ERROR\].*")) {
            stop(readText(log));        
        }    
    }

    return(0);
}

#' Create commandline arguments
#' 
#' @param argv a list object that contains commandline parameters for
#'    run a specific mothur command. commandline argument which it has 
#'    ``NULL`` value will be ignored in this function. 
#'
const mothurArgvs as function(command, argv) {
    for(name in names(argv)) {
        if (is.null(argv[[name]])) {
            warning(`missing value of '${name}' in '${command}'?`);
        }
    }

    argv 
    |> names()
    |> which(key -> !is.null(argv[[key]]))
    |> sapply(key -> `${key}=${argv[[key]]}`)
    |> paste(deli = ", ")
    ;
}

#' Parse error message
#' 
#' @param stdout the standard output dumpping text of the mothur
#'    program.
#' 
#' @return a error message pair list object.
#'  
const mothurErrorParser as function(stdout) {

}