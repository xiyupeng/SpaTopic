#' Calculate lift ranked words and values for a topic
#' 
#' @param beta A numeric matrix of dimension (topics x words) representing the probability
#' distribution of words within each topic. Each row should sum to 1. Beta must be on
#' the probability scale (not log scale).
#' @param vocab a character vector of vocabulary terms corresponding to the columns of beta.
#' @param wordcounts a numeric vector giving the total count of each word across the entire dataset.
#' @param topic the topic index that we want to calculate, the default is 1.
#' @param top_n the number of top words to return, the default is to return all words.
#' 
#' 
#' @return a data frame with ranks, words, and lift values of the words
#'
#' @export
lift_topic <- function(beta, vocab, wordcounts, topic = 1, top_n = NULL) {
  
  # overall word frequency
  word_freq <- wordcounts / sum(wordcounts)
  
  # lift
  lift_mat <- sweep(beta, 2, word_freq, "/")
  
  df <- data.frame(
    word = vocab,
    lift = lift_mat[topic, ]
  )
  
  # sort
  df <- df[order(df$lift, decreasing = TRUE), ]
  
  # limit output
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }
  
  rownames(df) <- NULL
  return(df)
}