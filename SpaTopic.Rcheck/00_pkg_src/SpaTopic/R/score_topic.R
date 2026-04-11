#' Calculate score ranked words and values for a topic
#' 
#' @param beta A numeric matrix of dimension (topics x words) representing the probability
#' distribution of words within each topic. Each row should sum to 1. Beta must be on
#' the probability scale (not log scale).
#' @param vocab a character vector of vocabulary terms corresponding to the columns of beta.
#' @param topic the topic index that we want to calculate, the default is 1.
#' @param top_n the number of top words to return, the default is to return all words.
#' 
#' 
#' @return a data frame with ranks, words, and score values of the words
#'
#' @export
score_topic <- function(beta, vocab, topic = 1, top_n = NULL) {
  
  beta <- pmax(beta, 1e-12) # avoid log(0)
  
  # score
  log_beta <- log(beta)
  avg_log_beta <- colMeans(log_beta)
  score_mat <- beta * (
    log_beta - matrix(avg_log_beta, nrow = nrow(beta), ncol = ncol(beta), byrow = TRUE)
  )
  
  df <- data.frame(
    word = vocab,
    score = score_mat[topic, ]
  )
  
  # sort
  df <- df[order(df$score, decreasing = TRUE), ]
  
  # limit output
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }
  
  rownames(df) <- NULL
  return(df)
}