replaceTabs.inner <-  function( text, width=8 ) # CH: For ap_textplot.character, copied from online, belongs to gplots (original texplot package), but seems to not be exported.
{
  spaces <- "        "

  if(nchar(text)<1) return(text)

  text.split <- strsplit(text,"\t")[[1]]
  if(length(text.split)==1)
    return(text)
  else
  {
    nSpaceAdd <- 8 - nchar(text.split) %% 8
    nSpaceAdd[length(nSpaceAdd)] <- 0
    nSpaceAdd[nSpaceAdd==8] <- 0

    retval <- ""
    for(i in 1:length(text.split))
    {
      tmp.text <- chartr("\t"," ", text.split[i]) # one space here
      retval <- paste(retval, tmp.text, substr(spaces,0,nSpaceAdd[i]-1 ), sep='' ) # rest here
    }
    return(retval)
  }
}

replaceTabs <- function(text, width=8) # CH: For ap_textplot.character, copied from online, belongs to gplots (original texplot package), but seems to not be exported.
{
  text <- as.character(text)
  retval <- sapply(text, replaceTabs.inner)
  names(retval) <- names(text)
  retval
}
