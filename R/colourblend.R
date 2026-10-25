#' Opaque colour blending
#'
#' Produce the opaque colour that looks like a semi-transparent colour drawn over a background, by mixing the two colours. This is useful for graphics devices that do not support semi-transparency, such as \code{postscript} (EPS), where semi-transparent elements are not drawn.
#'
#' @param col colours to blend. Any form accepted by \code{\link[grDevices]{col2rgb}}. Their own alpha channel is ignored.
#' @param alpha opacity of \code{col}, between 0 (only the background) and 1 (only \code{col}). \code{col} and \code{alpha} are recycled to a common length.
#' @param bg background colour. If \code{NULL} (default), the background colour of the current device, \code{par("bg")}, is used if a device is open, otherwise white. A transparent background is treated as white.
#'
#' @return A character vector of opaque colours in hexadecimal format.
#'
#' @details
#' Each colour channel is mixed as \code{alpha * col + (1 - alpha) * bg}, which is how a semi-transparent colour is drawn over an opaque background. The result is identical to using \code{adjustcolor(col, alpha.f=alpha)} on that background, but overlapping elements do not add up: the colour drawn last covers the ones below it. If the plot is placed on a different background afterwards, the blended colours keep the tint of \code{bg}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
#'
#' @seealso \code{\link[grDevices]{adjustcolor}}
#'
#' @examples
#' colour.blend("red", 0.25)
#' colour.blend("steelblue", c(0.2, 0.5, 0.8))
#' colour.blend(c("red", "blue"), 0.5, bg="black")
#' plot(NA, xlim=c(0,1), ylim=c(0,1))
#' rect(0.1, 0.1, 0.9, 0.9, col=colour.blend("steelblue", 0.3))
#'
#' @export colour.blend

colour.blend <- function(col, alpha, bg=NULL){

  if (length(col) == 0 || length(alpha) == 0){
    return(character(0))
  }
  if (any(is.na(alpha)) || any(alpha < 0) || any(alpha > 1)){
    stop("alpha must be between 0 and 1.")
  }

  # Background: current device if one is open (to avoid opening a new one), otherwise white
  if (is.null(bg)){
    bg <- if (grDevices::dev.cur() > 1) graphics::par("bg") else "white"
  }
  if ((length(bg) != 1) || is.na(bg) || (grDevices::col2rgb(bg, alpha=TRUE)[4] == 0)){
    if (length(bg) != 1){
      warning("Only the first background colour is used.")
      return(colour.blend(col, alpha, bg[1]))
    }
    bg <- "white"
  }

  # Mix channels
  n <- max(length(col), length(alpha))
  col <- rep_len(col, n)
  alpha <- rep_len(alpha, n)
  a <- grDevices::col2rgb(col)/255
  b <- grDevices::col2rgb(bg)[,1]/255
  w <- matrix(alpha, nrow=3, ncol=n, byrow=TRUE)
  m <- w*a + (1-w)*b
  grDevices::rgb(m[1,], m[2,], m[3,])

}
