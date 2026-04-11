#' Calculate FREX ranked words and values for a topic
#' 
#' @param beta A numeric matrix of dimension (topics x words) representing the probability
#' distribution of words within each topic. Each row should sum to 1. Beta must be on
#' the probability scale (not log scale).
#' @param vocab a character vector of vocabulary terms corresponding to the columns of beta.
#' @param topic the topic index that we want to calculate, the default is 1.
#' @param top_n the number of top words to return, the default is to return all words.
#' @param frex_weight the weight between 0 and 1 controlling the balance between
#' frequency and exclusivity in the FREX metric. Weight closer to 1 is favoring exclusivity
#' and closer to 0 is favoring frequency, we set the default as 0.5.
#' 
#' @return a data frame with ranks, words, and FREX values of the words
#'
#' @export
frex_topic <- function(beta, vocab, topic = 1, top_n = NULL, frex_weight = 0.5) {
  
  # frequency rank within each topic
  freq_rank <- t(apply(beta, 1, function(x) rank(x, ties.method = "average") / length(x)))
  
  # exclusivity
  ex <- sweep(beta, 2, colSums(beta), "/")
  ex_rank <- t(apply(ex, 1, function(x) rank(x, ties.method = "average") / length(x)))
  
  # FREX
  w <- frex_weight
  frex_mat<- 1 / (w / freq_rank + (1 - w) / ex_rank)
  
  df <- data.frame(
    word = vocab,
    frex = frex_mat[topic, ]
  )
  
  # sort
  df <- df[order(df$frex, decreasing = TRUE), ]
  
  # limit output number
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }
  
  rownames(df) <- NULL
  return(df)
}