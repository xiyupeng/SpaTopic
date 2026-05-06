#' Label topics function
#'
#' This function generates topic labels using four metrics:
#' highest probability, FREX, lift, and score. For each topic, it returns
#' the top n vocabulary terms according to each metric.
#'
#' @param beta A numeric matrix of dimension (topics x words) representing the probability
#' distribution of words within each topic. Each row should sum to 1. Beta must be on
#' the probability scale (not log scale).
#' @param vocab a character vector of vocabulary terms corresponding to the columns of beta.
#' @param wordcounts a numeric vector giving the total count of each word across the entire dataset.
#' @param n the number of top words to return for each topic, the default value is 8.
#' @param frex_weight the weight between 0 and 1 controlling the balance between
#' frequency and exclusivity in the FREX metric. Weight closer to 1 is favoring exclusivity
#' and closer to 0 is favoring frequency, we set the default as 0.5.
#'
#' @return a list of top n vocabulary terms for each topic, ranked according to four metrics:
#' highest probability, FREX, lift, and score.
#' 
#' @details
#' Highest Probability:
#'    For each topic, words are ranked by their probability within that topic.
#'    The top n words with the largest probabilities are selected.
#' FREX:
#'    FREX is calculated by combining frequency and exclusivity for each word in each topic.
#'    Frequency is the word probabilities ranked and scaled to values between 0 and 1.
#'    Each word’s probability is divided by its total probability to calculate how exclusive
#'    the word is to each topic. Then the exclusivity values are ranked within each topic and
#'    scaled to values between 0 and 1. The FREX score is the weighted harmonic mean of frequency
#'    rank and exclusivity rank, according to this formula frex<- 1 / (w / freq_rank + (1 - w) / ex_rank).
#' Lift:
#'    We first calculate the overall frequency of each word by dividing its total count
#'    by the total count of all words in the dataset. Then each word’s probability is divided by
#'    its overall frequency.
#' Score:
#'    The score is computed by first taking the logarithm of the topic-word probabilities.
#'    Then calculate the average log probability across all topics for each word to represent
#'    its overall baseline level. For each topic and word, compute the difference between
#'    its log probability in that topic and its average log probability, and multiply by
#'    beta to get the final score.
#'    
#' @export
label_topics <- function(beta, vocab, wordcounts, n = 8, frex_weight = 0.5) {
  K <- nrow(beta) # number of topics
  
  # helper function: get top n vocabs from rows
  top_n <- function(row, vocab, n) {
    vocab[order(row, decreasing = TRUE)][1:n]
  }
  
  # Highest Probability
  prob_labels <- apply(beta, 1, top_n, vocab, n)
  
  # FREX
  freq_rank <- t(apply(beta, 1, function(x) rank(x, ties.method = "average") / length(x)))
  ex <- sweep(beta, 2, colSums(beta), "/")
  ex_rank <- t(apply(ex, 1, function(x) rank(x, ties.method = "average") / length(x)))
  w <- frex_weight
  frex_mat<- 1 / (w / freq_rank + (1 - w) / ex_rank)
  
  frex_labels <- apply(frex_mat, 1, top_n, vocab, n)
  
  # Lift
  word_freq <- wordcounts/sum(wordcounts) # use freq instead of wordcounts, output should be the same
  lift_mat <- sweep(beta, 2, word_freq, "/")
  lift_labels <- apply(lift_mat, 1, top_n, vocab, n)
  
  # Score
  beta <- pmax(beta, 1e-12) # avoid log(0)
  log_beta <- log(beta)
  avg_log_beta <- colMeans(log_beta)
  score_mat <- beta * (log_beta - matrix(avg_log_beta, nrow = nrow(beta), ncol = ncol(beta), byrow = TRUE))
  score_labels <- apply(score_mat, 1, top_n, vocab, n)
  
  # results
  list(
    prob = prob_labels,
    frex = frex_labels,
    lift = lift_labels,
    score = score_labels
  )
}