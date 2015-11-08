# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

mat2wd = function(file_name) {
  f = file(file_name,open="r")
  lines = readLines(f)
  info_list = as.numeric(unlist(strsplit(lines[1], " ")))
  num_doc = info_list[1]
  num_vocab = info_list[2]
  nnz = info_list[3]
  WD = matrix(0, nrow = num_vocab, ncol = num_doc)
  
  for (i in 2:length(lines)) {
    l = as.numeric(unlist(strsplit(lines[i], " ")))
    for (j in seq(1, length(l), 2)) {
      word_index = l[j]
      word_occurence = l[j+1]
      WD[word_index, i-1] = word_occurence
    }
  }
  close(f)
  
  return(WD)
}

read_vocab = function(file_name) {
  f = file(file_name,open="r")
  vocab_list = readLines(f)
  close(f)
  return(vocab_list)
}

sparse_mat_count = function(WD) {
  nnz_indices = as.matrix(which(WD != 0, arr.ind = T))
  nnz = sum(WD)
  WS = rep(0, nnz)
  DS = rep(0, nnz)
  
  start_index = 1
  for(i in 1:dim(nnz_indices)[1]) {
    word_occurence = WD[nnz_indices[i,1], nnz_indices[i,2]] 
    WS[start_index:(start_index+word_occurence-1)] = nnz_indices[i, 1]
    DS[start_index:(start_index+word_occurence-1)] = nnz_indices[i, 2]
    start_index = start_index + word_occurence
  }
  result = list(ws=WS, ds=DS, nnz=nnz)
  return(result)
}

GibbsSamplerLDA = function(mat_txt, vocab_txt, k, alpha, beta, N){
  
    WD = mat2wd(mat_txt)
    g(W, D) %=% dim(WD)
    g(WS, DS, nnz) %=% sparse_mat_count(WD)
    WO = read_vocab(vocab_txt)
    V = length(WO)

    # initialize Z; store Z using names() and matrix
    #names(z_states) = 1:nnz
    z_states = matrix(data=-1, nrow=nnz, ncol=N)

    # initialize random topic assignment
    for (i in 1:nnz){
        z_states[i,1] = sample(1:k, size=1)
    }
    
    # update each word sequentially
    for (iter_ in 1:N){
        cat('Iteration: ', iter_)
        print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
      
        for (i in 1:nnz){
          # find the current word  wi
          wi = WS[i]
          # find all wi's in WS by indices
          wi_indices = array(which(WS %in% wi, arr.ind = TRUE))
          wi_topics_assignments = rep(-1, length(wi_indices))
          # find their up-to-date topic assignment
          for (index in 1:length(wi_indices)){
              wi_index = wi_indices[index]
              newest_index = tail(which(z_states[wi_index, ]!=-1, T), 1)
              wi_topics_assignments[index] = z_states[wi_index, newest_index]
          }
          
          topic_assignments = rep(-1, length(nnz))
          for (index in 1:nnz){
              # find up-to-date topic assignment for all nnz tokens
              newest_index = tail(which(z_states[index, ]!=-1, T), 1)
              topic_assignments[index] = z_states[index, newest_index]
          }
          
          # find the current doc di
          di = DS[i]
          # find all words of di in DS by indices
          w_di_indices = array(which(DS %in% di, arr.ind = TRUE))
          w_di_topics_assignments = rep(-1, length(w_di_indices))
          for (index in 1:length(w_di_indices)){
              di_index = w_di_indices[index]
              # find their up-to-date topic assignment
              newest_index = tail(which(z_states[di_index, ]!=-1, T), 1)
              w_di_topics_assignments[index] = z_states[di_index, newest_index]
          }
          
          zi_topic_prob = list()
          for (j in 1:k){
             newest_index = tail(which(z_states[i, ]!=-1, T), 1)
             # number of wi assigned to topic j excluding the current one
             n_wi_i_j = sum(wi_topics_assignments==j) - (z_states[i, newest_index] == j)
             # number of all words assigned to topic j excluding the current one
             n_wi_j = sum(topic_assignments==j) - (z_states[i, newest_index] == j)
             # number of words in di assigned to topic j excluding the current one
             n_di_i_j = sum(w_di_topics_assignments==j) - (z_states[i, newest_index] == j)
             
             # unnormalized probability
             zi_topic_prob = c(zi_topic_prob, (n_wi_i_j+beta)/(n_wi_j+W*beta) * (n_di_i_j+alpha))
          }
          
          # normalize probability
          zi_topic_normed_prob = unlist(zi_topic_prob)/sum(unlist(zi_topic_prob))
          # sample
          zi_sample = sample(1:k, size=1, prob=zi_topic_normed_prob)
          # update z
          z_states[i, iter_] = zi_sample
          
        }
    }
    
    # discard the first half as burn in
    approximate_states = z_states[,as.integer(N/2):N]
    
    # get topic distributions
    topic_distribution = c()
    
    # get topic-word distribution: matrix with shape(k, V)
    phi = matrix(data=-1, nrow = k, ncol = V)
    
    for (t in 1:k){
      # total times of words assigned to Topit t
      n_t = length(which(approximate_states==t, T))/2
      # topic proportion
      topic_distribution = c(topic_distribution, n_t/length(approximate_states))
      
      # total number of times w assigned to Topic t
      for (w in 1:V){
        indices = which(WS==w, T)
        n_w_t = length(which(approximate_states[indices, ]==t, T))/2
        phi[t,w] = (n_w_t + beta)/(n_t+W*beta)
      }
    }
    
    # get doc-topic distribution: matrix with shape(D, k)
    theta = matrix(data=-1, nrow = D, ncol = k)
    
    for (d in 1:D){
      # total times of words in Document d
      n_d = length(which(DS == d, T))
      
      # total number of times d assigned to Topic t
      for (t in 1:k){
        indices = which(DS==d, T)
        n_d_t = length(which(approximate_states[indices, ]==t, T))/2
        theta[d,t] = (n_d_t + alpha)/(n_d + k*alpha)
      }
    }
    
    print(topic_distribution)
    
    for (t in 1:k) {
      top_indices = sort(phi[t,], decreasing=T, index.return=T)$ix
      top_words = WO[top_indices[1:10]]
      writeLines(cat(sprintf("Topic: %i", t), top_words))
    }
    
    
    for (d in 1:D) {
      top_indices = sort(theta[d,], decreasing=T, index.return=T)$ix
      writeLines(cat(sprintf("Doc %i:", d), top_indices))
    }
}
  
LDA = function() {
  args <- commandArgs(trailingOnly = TRUE)
  GibbsSamplerLDA(args[1], args[2], as.numeric(args[3]), as.numeric(args[4]), as.numeric(args[5]), as.numeric(args[6]))
}

LDA()
