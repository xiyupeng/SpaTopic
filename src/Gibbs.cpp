//#include <Rcpp.h>


#include <cstdlib>
#include <Rmath.h>
#include <algorithm>
#include <fstream>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;


size_t discrete_sample(double * P, size_t K)
{
  size_t i = 0;
  for (i = 1; i < K; i++) {
    P[i] += P[i - 1];
  }
  double alpha = Rf_runif(0.0, P[K - 1]);
  double *target = std::upper_bound(P, P + K, alpha);
  return (size_t)(target - P);
}

// [[Rcpp::export]]
IntegerVector table_1d_fast(IntegerVector input, int n)
{
  IntegerVector res(n);
  
  int n_pairs = input.size();
  for (int i = 0; i < n_pairs; i++) {
    ++res(input(i));
  }
  return res;
}

// [[Rcpp::export]]
IntegerMatrix table_2d_fast(IntegerVector lhs, IntegerVector rhs, int n_lhs, int n_rhs)
{
  IntegerMatrix res(n_lhs, n_rhs);
  int n_pairs = lhs.size();
  for (int i = 0; i < n_pairs; i++) {
    ++res(lhs(i), rhs(i));
  }
  return res;
}

arma::uvec gen_uvec(IntegerVector src)
{
  arma::uvec dst(src.size());
  int i = 0;
  for (i = 0; i < src.size(); i++) {
    dst(i) = (arma::uword)src(i);
  }
  return dst;
}

// [[Rcpp::export]]
double compute_loglike(arma::mat m_theta, arma::mat m_beta, IntegerMatrix docs, IntegerMatrix neighbors, NumericMatrix Kernel, size_t M, 
size_t n_words, size_t K, double beta = .05, double alpha = .01, double sigma = 50)
{

  size_t n_pairs = docs.nrow();
  size_t n_neighbor = neighbors.ncol();
  
  size_t j, k, nc;
  //double sum;
  double perp = 0.;

  //arma::mat m_theta(theta, M, K, false);
  //arma::mat m_beta(beta, n_words, K, false);
  //arma::mat m_eta(n_pairs,n_neighbor);
  arma::mat loglike(n_pairs, K);
  arma::vec ext(M);
 

  for (j = 0; j < n_pairs; j++) {

    // Initialize with zero
    //sum = 0.;
    ext.zeros(M);

    // P(x_gi^c, Dgi = d | x^d_d)
    //for(nc = 0; nc < n_neighbor; nc++){
    //  m_eta(j,nc) = exp(-dists(j,nc)/sigma);
    //  sum += m_eta(j,nc);
    //}

    // normalizing to one and reformat
    for (nc = 0; nc < n_neighbor; nc++){
      //m_eta(j,nc) /= sum;
      ext(neighbors(j,nc)) = Kernel(j,nc);
    }

    // \sum_d P(Zgi = k, Dgi = d, x_gi^c | x^d_d)
    for (k = 0; k < K; k++){
      loglike(j,k) = accu(ext % m_theta.col(k));
    }

    // \sum_k \sum_d P(Zgi = k, Dgi = d, x_gi^c, C_gi | x^d_d)
    loglike.row(j) = m_beta.row(docs(j,1)) % loglike.row(j);
    perp += log(accu(loglike.row(j)));

  }

  //perp = exp(-perp/n_pairs);

  return(perp);
}


// [[Rcpp::export]]
List gibbs_lda_c(IntegerMatrix docs, IntegerMatrix Ndk, IntegerMatrix Nwk, IntegerVector Nk, IntegerVector Z,IntegerVector D,
 IntegerVector doc_list, IntegerVector word_list, size_t K, size_t niter = 20, double beta = .05, double alpha = .01)
{
  size_t n_pairs = docs.nrow();
  size_t i, j, k;
  double *P;
  
  size_t n_words = Nwk.nrow();
  
 // Rprintf("IN the function now! 1\n");
  
  arma::icube doc_trace(doc_list.size(), K, niter);
  arma::icube word_trace(word_list.size(), K, niter);
  arma::imat Nk_trace(K, niter);
  
 // Rprintf("IN the function now! 2\n");

  arma::icolvec m_Z(Z.begin(), Z.size(), false);
  arma::imat m_Ndk(Ndk.begin(), Ndk.nrow(), Ndk.ncol(), false);
  arma::imat m_Nwk(Nwk.begin(), Nwk.nrow(), Nwk.ncol(), false);
  arma::icolvec m_Nk(Nk.begin(), Nk.size(), false);

   // Rprintf("IN the function now! 3\n");
  
  arma::uvec m_doc_list = gen_uvec(doc_list);
  arma::uvec m_word_list = gen_uvec(word_list);

   // Rprintf("IN the function now! 4\n");
  
  GetRNGstate();
  
  P = new double[K];
  
  for (i = 0; i < niter; i++) {
    for (j = 0; j < n_pairs; j++) {
      
    
    // Rprintf("niter %ld; j: %ld of %ld \n",i,j,n_pairs);

      size_t doc_id = D(j);
      size_t word_id = docs(j, 1);
      //size_t cnt = docs(j, 2);
      size_t topic = m_Z(j);
    // Rprintf("%ld  %ld  %ld  %ld \n",m_Z(1),m_Z(2),m_Z(3),m_Z(4));
      
    // Rprintf("doc_id %ld; word_id: %ld; topic: %ld \n",doc_id,word_id,topic);

      --m_Ndk(doc_id, topic);
      --m_Nwk(word_id, topic);
      --m_Nk(topic);
     // --m_Nd(doc_id);
      
      for (k = 0; k < K; k++) {
        P[k] = (m_Ndk(doc_id, k) + alpha) * (m_Nwk(word_id, k) + beta) / (m_Nk(k) + beta * n_words);
      }
     // Rprintf("%f  %f  %f \n",P[0],P[1],P[2]);

      topic = discrete_sample(P, K);
      
     //  Rprintf("doc_id %ld; word_id: %ld; topic: %ld \n",doc_id,word_id,topic);

      m_Z(j) = topic;
      ++m_Ndk(doc_id, topic);
      ++m_Nwk(word_id, topic);
      ++m_Nk(topic);
      // ++m_Nd(doc_id);
    }

    Nk_trace.col(i) = m_Nk;

    //Rprintf("1");

    doc_trace.slice(i) = m_Ndk.rows(m_doc_list);

    //Rprintf("2");

    word_trace.slice(i) = m_Nwk.rows(m_word_list);

    //Rprintf("3");
  }
  
  PutRNGstate();
  
  delete [] P;
  
  return Rcpp::List::create(Rcpp::Named("doc.trace") = doc_trace,
                            Rcpp::Named("word.trace") = word_trace,
                            Rcpp::Named("nk.trace") = Nk_trace
  );
}

int gibbs_spaTopic_dev(IntegerMatrix docs, arma::imat& m_Ndk, arma::imat& m_Nwk, arma::icolvec& m_Nk, arma::icolvec& m_Nd,
 arma::icolvec& m_Z, arma::icolvec& m_D, IntegerMatrix neighbors,  NumericMatrix Kernel, size_t M, size_t n_words, size_t K, size_t niter = 20, 
  double beta = .05, double alpha = .01, double sigma = 50)
{
  size_t n_pairs = docs.nrow();
  //size_t M = m_Ndk.nrow();
  //size_t n_words = m_Nwk.nrow();
  size_t n_neighbor = neighbors.ncol();
  
  size_t i, j, k, nc;
  double *Pz;
  double *Pd;

  
  //arma::icolvec m_Z(Z, n_pairs, false);
  //arma::icolvec m_D(D, n_pairs, false);
  //arma::imat m_Ndk(Ndk, M, K, false);
  //arma::imat m_Nwk(Nwk, n_words, K, false);
  //arma::icolvec m_Nk(Nk, K, false);
  //arma::icolvec m_Nd(Nd, M, false);
  //arma::icolvec m_Z = &Z;
  //arma::icolvec m_D = &D;
  //arma::imat m_Ndk = &Ndk;
  //arma::imat m_Nwk = &m_Nwk;
  //arma::icolvec m_Nk = &m_Nk;
  //arma::icolvec m_Nd = &m_Nd;
 
  
  Pz = new double[K];
  Pd = new double[n_neighbor];

  double beta_nwords = beta * n_words;
  double alpha_K = alpha * K;
  // double sigma_2 = sigma * sigma;
  
  for (i = 0; i < niter; i++) {

    Rcpp::checkUserInterrupt();
    
    for (j = 0; j < n_pairs; j++) {
      
    
    // Rprintf("niter %ld; j: %ld of %ld \n",i,j,n_pairs);

    //size_t image_id = docs(j, 0);  // make the image id in the first column
      size_t word_id = docs(j, 1);   // the word id in the second column
      size_t doc_id = m_D(j);       // D is changing during iterations 
      size_t topic = m_Z(j);        // Z is changing during iterations
    // Rprintf("%ld  %ld  %ld  %ld \n",m_Z(1),m_Z(2),m_Z(3),m_Z(4));
      
   // Rprintf("doc_id %ld; word_id: %ld; topic: %ld \n",doc_id,word_id,topic);

      --m_Ndk(doc_id, topic);
      --m_Nwk(word_id, topic);
      --m_Nk(topic);
      --m_Nd(doc_id);
      
      // Update topic
      for (k = 0; k < K; k++) {
        Pz[k] = (m_Ndk(doc_id, k) + alpha) * (m_Nwk(word_id, k) + beta) / (m_Nk(k) + beta_nwords);
      }
      // Rprintf("%f  %f  %f \n",Pz[0],Pz[1],Pz[2]);

      topic = discrete_sample(Pz, K); // need check before adding to the vector
      if(topic >  K)
        Rcpp::stop("Assigned topic id exceed number of existing topic");

      m_Z(j) = topic;
      // Rprintf("doc_id %ld; word_id: %ld; topic: %ld \n",doc_id,word_id,topic);

      // update document
      // prob of cells to the nearest neighbor centers
      double max_prob = -1.;
      for (nc = 0; nc < n_neighbor; nc++){
        size_t doc = neighbors(j,nc);  //docs
        //double d = dists(j,nc);
       // Rprintf("%ld \n",m_Nd(doc));
        double prob = (m_Ndk(doc, topic) + alpha)/(m_Nd(doc)+ alpha_K) * Kernel(j,nc);
      //  Rprintf("j %ld; nc %ld; doc: %ld; d: %f prob: %f \n",j, nc, doc, d, prob);
        Pd[nc] = prob;
        if(prob > max_prob) max_prob = prob;
      }

      // avoid tiny small numbers
      double scale_factor = 1/(max_prob);
      for (nc = 0; nc < n_neighbor; nc++)
        Pd[nc] = Pd[nc] * scale_factor;
        
        
      // Rprintf("%f  %f  %f %f  %f\n",Pd[0],Pd[1],Pd[2],Pd[3],Pd[4]);
      nc = discrete_sample(Pd, n_neighbor); // need check before add to the vector
      if(nc >  n_neighbor)
        Rcpp::stop("Assigned neighbor id exceed number of identified neighbors");

      doc_id = neighbors(j,nc);

      if(doc_id >  M)
        Rcpp::stop("Assigned region id exceed number of existing regions");
      m_D(j) = doc_id;

      
      ++m_Ndk(doc_id, topic);
      ++m_Nwk(word_id, topic);
      ++m_Nk(topic);
      ++m_Nd(doc_id);
    }
  }
 
  delete [] Pz;
  delete [] Pd;

  return 0;
}

// [[Rcpp::export]]
Rcpp::List gibbs_sampler_c(IntegerMatrix docs, IntegerMatrix Ndk, IntegerMatrix Nwk, IntegerVector Nk, IntegerVector Nd, 
  IntegerVector Z, IntegerVector D, NumericMatrix dists,  IntegerMatrix neighbors, IntegerVector doc_list, IntegerVector word_list, 
  size_t K, double beta = .05, double alpha = .01, double sigma = 50, size_t thin = 20, size_t burnin = 1000, size_t niter = 200, 
  bool trace = false,bool display_progress =true){

    size_t i,j,k,nc;
    size_t word_id,doc_id;
    size_t M = Ndk.nrow();
    size_t n_words = Nwk.nrow();
    size_t n_pairs = docs.nrow();
    size_t n_neighbor = neighbors.ncol();
    double sum;
    double *Pz;

    Progress p(niter, display_progress);
  
    arma::icube doc_trace;
    arma::icube word_trace;
   
    arma::mat Z_prob(n_pairs,K);
    arma::imat Z_trace; 
    arma::mat Theta;
    arma::mat Beta;
    
    
    NumericMatrix Kernel(neighbors.nrow(),neighbors.ncol());
    // consider the memory issue, only save the total count, other than the occurrence at each iteration
    Z_trace.zeros(Z.size(),K);
    Theta.zeros(Ndk.nrow(),Ndk.ncol());
    Beta.zeros(Nwk.nrow(),Nwk.ncol());
   
    

    if(trace){
      doc_trace.zeros(doc_list.size(), K, niter);
      word_trace.zeros(word_list.size(), K, niter);
    }

    arma::icolvec m_Z(Z.begin(), Z.size(), false);
    arma::icolvec m_D(D.begin(), D.size(), false);
    arma::imat m_Ndk(Ndk.begin(), Ndk.nrow(), Ndk.ncol(), false);
    arma::imat m_Nwk(Nwk.begin(), Nwk.nrow(), Nwk.ncol(), false);
    arma::icolvec m_Nk(Nk.begin(), Nk.size(), false);
    arma::icolvec m_Nd(Nd.begin(), Nd.size(), false);

    arma::uvec m_doc_list = gen_uvec(doc_list);
    arma::uvec m_word_list = gen_uvec(word_list);
    arma::vec loglike(niter+1);

    GetRNGstate();

    Pz = new double[K];

  // precompute the kernel function
    for(j = 0; j < n_pairs; j++){
      sum = 0.;
      for (nc = 0; nc < n_neighbor; nc++){
          double d = dists(j,nc);
          double prob = exp(-d /sigma);
          Kernel(j,nc) = prob;
          sum += prob;
      }
      // normalize 
      for (nc = 0; nc < n_neighbor; nc++)
        Kernel(j,nc)/= sum;
    }

    
    // burn-in
    gibbs_spaTopic_dev(docs, m_Ndk, m_Nwk, m_Nk, m_Nd, m_Z, m_D, 
                        neighbors, Kernel, M, n_words, K, burnin, beta, alpha, sigma);
    

    // Gibbs sampling after burn-in
    for (i = 0; i < niter; i++) {

      // Calculate current Beta and Theta
        if(trace){
          for (k = 0; k < K; k++) {
            for (word_id = 0; word_id < n_words; word_id++)
              Beta(word_id,k) = (m_Nwk(word_id, k) + beta) / (m_Nk(k) + beta*n_words);
            for (doc_id = 0; doc_id < M; doc_id++)  
              Theta(doc_id,k) = (m_Ndk(doc_id, k) + alpha)/(m_Nd(doc_id)+ alpha*K); 
          }
          loglike(i) = compute_loglike(Theta, Beta, docs, neighbors, Kernel, M, n_words, K, beta, alpha, sigma);  
        }

        //Rprintf("niter %ld; \n",i);
        p.increment(); // update progress
        gibbs_spaTopic_dev(docs, m_Ndk, m_Nwk, m_Nk, m_Nd, m_Z, m_D,
                          neighbors, Kernel, M, n_words, K, thin, beta, alpha, sigma);


        for(j = 0; j < n_pairs; j++){
          size_t topic_id = m_Z(j); 
          ++Z_trace(j,topic_id);
        }
            
        if(trace){
          doc_trace.slice(i) = m_Ndk.rows(m_doc_list);
          word_trace.slice(i) = m_Nwk.rows(m_word_list);
        }

    }

    // Calculate Beta and Theta from the final sample
      for (k = 0; k < K; k++) {
        for (word_id = 0; word_id < n_words; word_id++)
          Beta(word_id,k) = (m_Nwk(word_id, k) + beta) / (m_Nk(k) + beta*n_words);
        for (doc_id = 0; doc_id < M; doc_id++)  
          Theta(doc_id,k) = (m_Ndk(doc_id, k) + alpha)/(m_Nd(doc_id)+ alpha*K); 
      }

    // calculate the final model predictive log likelihood
    loglike(niter) = compute_loglike(Theta, Beta, docs, neighbors, Kernel, M, n_words, K, beta, alpha, sigma);

    // calculate the model perplexity and deviance
    double perplexity = exp(-loglike(niter)/n_pairs);
    double deviance = -2*loglike(niter);
    // [TODO] calculate the convergence criterior of MCMC chain.

    PutRNGstate();

    delete [] Pz;
    
    return Rcpp::List::create(Rcpp::Named("Perplexity") =  perplexity,
                              Rcpp::Named("Deviance") =  deviance,
                              Rcpp::Named("loglikelihood") =  loglike(niter),
                              Rcpp::Named("loglike.trace") = loglike,
                              Rcpp::Named("Beta")  = Beta,
                              Rcpp::Named("Theta") = Theta,
                              Rcpp::Named("Ndk")  = m_Ndk,
                              Rcpp::Named("Nwk") = m_Nwk,
                              Rcpp::Named("Z.trace") = Z_trace,
                              Rcpp::Named("doc.trace") = doc_trace,
                              Rcpp::Named("word.trace") = word_trace);                           
  
}

