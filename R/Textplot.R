#' Display text information in a graphics plot.
#'
#' This function displays text output in a graphics window.
#' There is an option to display the text using the largest font that will fit in the plotting region.
#' For matrixes, data.frames and vectors a specialized textplot function
#' is available that plots each of the cells individually in a way that is visually
#' appealing (maintains the table-like/grid-like structure of the data).
#' If present, row and column labels will be displayed in a bold font.
#'
#' This function was modified from the PerformanceAnalytics version by Peter Carl and Brian G. Peterson (brian@@braverock.com).
#' Assistance was provided by John Colby in this post:
#' \url{http://stackoverflow.com/questions/8523944/left-justify-a-column-using-textplot-gplots-or-performanceanalytics}
#'
#' Modifications:
#' \itemize{
#'             \item More descriptive variables names, more in-line comments,
#'             more detailed/clear parameter documentation, more concise code
#              \item This roxygen2 documentation
#'             \item Bug fix and easier to use/understand cex parameter.
#'             If cex is specified then what you specify is what you get.
#'             Previously a specified cex was adjusted by par(cex) by the text() function, which is confusing.
#'             \item This function sets par(cex = 1) and then passes the cex parameter
#'             (or computed cex in the auto-size case) to text(), so caller simply decides
#'             what cex is desired or omits cex to auto-size.
#'             \item Bug fix to overlapping text when left aligning matrix/data.frame
#'             (see StackOverflow post linked above)
#'             \item Eliminate vadj and hadj in favor of vCellAlign and hCellAlign and renamed
#'             valign/halign to vPlotAlign and hPlotAlign.
#'             Caller decides "center", "left", "right", "top", "bottom" for plot and cells and it all works out.
#'             \item Using a small hack, textplot.character will horizontally align text within
#'             the plot region and within the plot itself.
#'             Previously the text was only aligned within the plot region so it appeared
#'             as if the alignment instruction was ignored.
#' }
#'
#' Color has NOT been tested.
#' Don't deal with the TAB character.
#'
#'
#' @param object            Object to be displayed.  If character of length > 1, then elements
#'     are concatenated with a newline delimiter.
#' @param hPlotAlign        Character with alignment of the entire object within the plotting region
#'     in the x direction, one of "center", "left", or "right".  For character \samp{object} with multiple lines of text
#'     this parameter is also used to align the text within the plot
#'     (i.e. "center" will center the object within the plotting region and also center 2
#'     or more lines of text within the plot itself)
#' @param vPlotAlign        Character with alignment of the entire object within the plotting region in the y direction,
#'     one of "center", "left", or "right".
#' @param cex               Numeric of length 1 with the text character size, see \link{par} for details.
#'     If missing, the code will attempt to use the largest value which allows the entire plot region to be filled.
#' @param cmar              Numeric with the space between columns, in fractions of the size of the letter 'M'.  Typically > 1
#' @param rmar              Numeric with the space between columns, in fractions of the size of the letter 'M'.  Typically > 1
#' @param showRowNames      Logical indicating whether row names will be displayed (for matrix, data.frame, or vector objects)
#' @param showColumnNames   Logical indicating whether column names will be displayed (for matrix, data.frame, or vector objects)
#' @param hCellAlign        Character with the horizontal location of elements within cell
#'     (for matrix, data.frame, or vector objects).  One of "center", "left", or "right".
#' @param vCellAlign        Character with the vertical location of elements within cell
#'     (for matrix, data.frame, or vector objects).  One of "center", "top", or "bottom".
#' @param mar               Figure margins, see the documentation for \link{par}
#' @param dataColor         Colors for data elements. If a single value is provided,
#'     all data elements will be the same color. If a matrix matching the dimensions of the data is provided,
#'    each data element will receive the specified color.
#' @param rowNameColor      Colors for row names. May be specified as a scalar or a vector of the same row length as object.
#' @param columnNameColor   Colors for column names. May be specified as a scalar or a vector of the same column length as object.
#' @param ...               Optional arguments passed to the text plotting command.
#'
#' @return    Numeric with the auto-sized cex value if cex is missing, otherwise the original cex value.
#'
#' @author Suraj Gupta (suraj at wingedfootcapital.com)
#' @export
Textplot = function( object , ... )
{
    UseMethod( "Textplot" )
}

#' @rdname   Textplot
#' @method   Textplot default
#' @S3method Textplot default
Textplot.default = function( object , ... )
{
    if ( is.matrix( object ) || ( is.vector( object ) && ( length( object ) > 1 ) ) )
    {
        return( Textplot.matrix( object , ... ) )
    }
    return( Textplot.character( as.character( object ) , ... ) )
}

#' @rdname   Textplot
#' @method   Textplot matrix
#' @S3method Textplot matrix
Textplot.matrix = function( object , hPlotAlign = c( "center" , "left" , "right" ) , vPlotAlign = c( "center" , "top" , "bottom" ),
                            cex , mar = c( 1 , 1 , 4 , 1 ) + 0.1 , cmar = 1 , rmar = 1 ,
                            showRowNames = TRUE , showColumnNames = TRUE , hCellAlign= c( "center" , "left" , "right" ) ,
                            vCellAlign = c( "center" , "top" , "bottom" ) , dataColor = par( "col" ) , rowNameColor = par( "col" ),
                            columnNameColor = par( "col" ) , ... )
{

    # convert vector to matrix
    if ( is.vector( object ) ) { object = t( as.matrix( object ) ) }
    else                       { object = as.matrix( object ) }

    # pull alignment parameters
    hPlotAlign = match.arg( hPlotAlign )
    vPlotAlign = match.arg( vPlotAlign )
    hCellAlign = match.arg( hCellAlign )
    vCellAlign = match.arg( vCellAlign )
    vadj = hadj = NULL
    if ( hCellAlign == "center" )    { hadj = .5 }
    else if ( hCellAlign == "left" ) { hadj = 0  }
    else                             { hadj = 1  }
    if ( vCellAlign == "center" )    { vadj = .5 }
    else if ( vCellAlign == "top" )  { vadj = 1  }
    else                             { vadj = 0  }

    # check dimensions of dataColor, rowNameColor, columnNameColor
    if( length( dataColor ) == 1 ) { dataColor = matrix( dataColor , nrow = nrow( object ) , ncol = ncol( object ) ) }
    else { if ( nrow( dataColor ) != nrow( object ) || ncol( dataColor ) != ncol( object ) ) {
      stop( "Dimensions of 'dataColor' do not match dimensions of 'object'." ) } }
    if( length( rowNameColor ) == 1 ) { rowNameColor = rep( rowNameColor , nrow( object ) ) }
    else { if ( length( rowNameColor ) != nrow( object ) ) {
      stop( "Length of 'rowNameColor' do not match number of rows in 'object'." ) } }
    if( length( columnNameColor ) == 1 ) { columnNameColor = rep( columnNameColor , ncol( object ) ) }
    else { if ( length( columnNameColor ) != ncol( object ) ) {
      stop( "Length of 'columnNameColor' do not match number of rows in 'object'." ) } }

    # save old par settings
    opar = par()[ c( "mar" , "xpd" , "cex" ) ]
    on.exit( par( opar ) )

    # set margins and force plotting to be clipped to plotting region
    # here we HAVE to set cex = 1 because text() multiplies the cex parmeter by par( cex )
    # to get the final cex which is terribly confusing.
    par( mar = mar , xpd = FALSE , cex = 1 )

    # start the plot and set up the coordinate system which we'll use for alignment
    plot.new()
    plot.window( xlim = c( 0 , 1 ) , ylim = c( 0 , 1 ) , log = "" , asp = NA )

    # add row/column names into the matrix itself as first row/column
    # copy the row/column colors into dataColor matrix
    if ( showRowNames )
    {
        if ( !is.null( rownames( object ) ) )
        {
            object = cbind( rownames( object ) , object )
            dataColor = cbind( rowNameColor , dataColor )
        }
    }
    if ( showColumnNames )
    {
        if ( !is.null( colnames( object ) ) )
        {
            object = rbind( colnames(object), object )
            dataColor = rbind( columnNameColor , dataColor )
        }
    }

    # if cex is not provided, then auto-size the text to fill the plot area
    if ( missing( cex ) )
    {
        cex = 1
        for ( i in 1 : 20 )
        {

            # get the sum of the widths of each column (take the widest text in each column and sum the widths)
            # then pad each column with cmar
            width = sum( apply( object, 2 , function( x ) max( strwidth( x , cex = cex ) ) ) ) +
              strwidth( "M" , cex = cex ) * cmar * ncol( object )

            # the height is simply the height of any character times the number of rows * the row margin
            height = strheight( 'M' , cex = cex ) * nrow( object ) * rmar
            oldcex = cex
            cex = cex / max( width , height )
            if (abs( oldcex - cex ) < 0.001 ) { break }
        }
    }

    # alignment
    width = sum( apply( object, 2 , function( x ) max( strwidth( x , cex = cex ) ) ) ) +
      strwidth( "M" , cex = cex ) * cmar * ncol( object )
    height = strheight( 'M' , cex = cex ) * nrow( object ) * rmar
    if ( hPlotAlign == "left" )        { xpos = 0 }
    else if ( hPlotAlign == "center" ) { xpos = 0 + ( 1 - width ) / 2 }
    else                               { xpos = 0 + ( 1 - width ) }
    if ( vPlotAlign == "top" )         { ypos = 1 }
    else if ( vPlotAlign == "center" ) { ypos = 1 - ( 1 - height ) / 2 }
    else                               { ypos = 1 - ( 1 - height ) }


    # iterate across elements, plotting them
    y = ypos
    colWidths = apply( object, 2 , function( x ) max( strwidth( x , cex = cex ) ) ) + strwidth( "M" ) * cmar
    rowHeight = strheight( "W" , cex = cex ) * rmar
    for ( i in 1 : ncol( object ) )
    {
        xpos = xpos + hadj * colWidths[ i ]
        for( j in 1 : nrow( object ) )
        {
            ypos = y - ( j - 1 ) * rowHeight
            if ( ( showRowNames && i == 1 ) || ( showColumnNames && j == 1 ) )
            {
                text( xpos , ypos , object[ j , i ] , adj = c( hadj , vadj ) , cex = cex , font = 2 ,
                      col = dataColor[ j , i ] , ... )
            }
            else
            {
                text( xpos , ypos , object[ j , i ] , adj = c( hadj , vadj ) , cex = cex , font = 1 ,
                      col = dataColor[ j , i ] , ... )
            }
        }
        xpos = xpos + ( 1 - hadj ) * colWidths[ i ]
    }

    # return resulting cex
    invisible( cex )
}

#' @rdname   Textplot
#' @method   Textplot data.frame
#' @S3method Textplot data.frame
Textplot.data.frame = function( object , ... )
{
    return( Textplot.matrix( object , ... ) )
}

#' @rdname   Textplot
#' @method   Textplot character
#' @S3method Textplot character
Textplot.character = function( object , hPlotAlign = c( "center" , "left" , "right" ) ,
                               vPlotAlign = c( "center" , "top" , "bottom" ) , cex , mar = c( 0 , 0 , 3 , 0 ) + 0.1 , ... )
{

    # if length > 1, then combine with new line
    object = paste( object , collapse = "\n" , sep = "" )

    # TODO: took out the code to handle TAB characters

    # pull hPlotAlign, vPlotAlign
    hPlotAlign = match.arg( hPlotAlign )
    vPlotAlign = match.arg( vPlotAlign )

    # start the plot
    plot.new()

    # save old par settings
    opar = par()[ c( "mar" , "xpd" , "cex" ) ]
    on.exit( par( opar ) )

    # set margins and force plotting to be clipped to plotting region
    # here we HAVE to set cex = 1 because text() multiplies the cex parmeter by par( cex )
    # to get the final cex which is terribly confusing.
    par( mar = mar , xpd = FALSE , cex = 1 )

    # set up the coordinate system which we'll use for alignment
    plot.window( xlim = c( 0 , 1 ) , ylim = c( 0 , 1 ) , log = "" , asp = NA )

    # if cex is missing then auto-size
    lastloop = FALSE
    if ( missing( cex ) )
    {
        cex = 1
        for ( i in 1 : 20 )
        {
            width = strwidth( object , cex = cex )
            height = strheight( object , cex = cex )
            oldcex = cex
            cex = cex / max( width , height )
            if ( abs( oldcex - cex ) < 0.001 ) { break }
        }
    }

    # alignment
    width = strwidth( object , cex = cex )
    height = strheight( object , cex = cex )
    if ( hPlotAlign == "left" )        { xpos = 0 }
    else if ( hPlotAlign == "center" ) { xpos = 0 + ( 1 - width ) / 2 }
    else                               { xpos = 0 + ( 1 - width ) }
    if ( vPlotAlign == "top" )         { ypos = 1 }
    else if ( vPlotAlign == "center" ) { ypos = 1 - ( 1 - height ) / 2 }
    else                               { ypos = 1 - ( 1 - height ) }

    # horizonal alignment doesn't work as-is for multi-line text.  that's because the entire text region is being aligned
    # and not the text itself. A bit of a hack, but we simply pad the shorter lines of text with spaces to achieve the
    # effect of left and right justification
    if ( ( hPlotAlign == "right" ) || ( hPlotAlign == "center" ) )
    {
        textByLine = unlist( lapply( object , function( x ) strsplit( x , "\n" ) ) )
        for ( i in 1 : length( textByLine ) )
        {
            oneLineText = textByLine[ i ]
            while ( strwidth( oneLineText , cex = cex ) <= width )
            {
                if ( hPlotAlign == "right" ) { oneLineText = paste( " " , oneLineText , sep = "" ) }
                else                     { oneLineText = paste( " " , oneLineText , " " , sep = "" ) }
            }
            textByLine[ i ] = oneLineText
        }
        object = paste( textByLine , collapse = "\n" )
    }

    # plot the text
    text( x = xpos , y = ypos , labels = object , adj = c( 0 , 1 ) ,  cex = cex , ... )

    # return resulting cex
    invisible( cex )
}
