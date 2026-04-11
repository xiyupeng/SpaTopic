#' Calculate highest probability words and values for a topic
#' 
#' @param beta A numeric matrix of dimension (topics x words) representing the probability
#' distribution of words within each topic. Each row should sum to 1. Beta must be on
#' the probability scale (not log scale).
#' @param vocab a character vector of vocabulary terms corresponding to the columns of beta.
#' @param topic the topic index that we want to calculate, the default is 1.
#' @param top_n the number of top words to return, the default is to return all words.
#'
#' @return a data frame with ranks, words, and the probabilities of the words
#'
#' @export
prob_topic <- function(beta, vocab, topic = 1, top_n = NULL) {
  df <- data.frame(
    word = vocab,
    prob = beta[topic, ]
  )
  
  # sort
  df <- df[order(df$prob, decreasing = TRUE), ]
  
  # limit output number
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }
  
  df$rank <- seq_len(nrow(df))
  df <- df[, c("rank", "word", "prob")]
  
  rownames(df) <- NULL
  return(df)
}