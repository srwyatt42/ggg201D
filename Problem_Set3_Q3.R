library(ggplot2)


#Time to first coalescent event
E_coalesce <- function(gene_copies){
  E_coal_data <- NULL
  for(k in gene_copies){
    t = 2/((k*(k-1))/2)
    E_coal_data <- append(E_coal_data, t)
  }
  return(E_coal_data)
}

#Time to MRCA
E_MRCA <- function(gene_copies){
  E_MRCA_data <- NULL
  for(n in gene_copies){
    k_values <- 2:n
    t_values <- NULL
    for(k in k_values){
      t = 2/((k*(k-1))/2)
      t_values <- append(t_values, k*t) 
    }
  E_MRCA_data <- append(E_MRCA_data, sum(t_values))
  }
  return(E_MRCA_data)
}

#Total Tree Length
E_TTL <- function(gene_copies){
  E_TTL_data <- NULL
  for(n in gene_copies){
    k_values <- 2:n
    t_values <- NULL
    for(k in k_values){  
    t = (2*k)/((k*(k-1))/n)
    t_values <- append(t_values, t)
    }
    E_TTL_data <- append(E_TTL_data, sum(t_values))
  }
  return(E_TTL_data)
}

E_TTL2 <- function(gene_copies){
  E_TTL2_data <- NULL
  for(n in gene_copies){
    i_values <- 1:(n-1)
    t_values <- NULL
    for(i in i_values){
      t = 1/i
      t_values <- append(t_values, t)
    }
    E_TTL2_data <- append(E_TTL2_data, 2*sum(t_values))
  }
  return(E_TTL2_data)
}

genecopies = 2:50
coalesce = E_coalesce(genecopies)
MRCA = E_MRCA(genecopies)
TTL = E_TTL(genecopies)
TTL2 = E_TTL2(genecopies)

Data <- data.frame(genecopies, coalesce, MRCA, TTL, TTL2)

#Graph
ggplot(Data, aes(genecopies, coalesce)) +
  geom_point(color = "purple") +
  labs(x = "Number of Gene Copies", y = "Expected Number of Generations in N units", title = "Expected Value of Time to Coalescence")

ggplot(Data, aes(genecopies, MRCA)) +
  geom_point(color = "purple") +
  labs(x = "Number of Gene Copies", y = "Expected Number of Generations in N units", title = "Expected Value of Time to MRCA")

ggplot(Data, aes(genecopies, TTL)) +
  geom_point(color = "purple") +
  labs(x = "Number of Gene Copies", y = "Expected Number of Generations in N units", title = "Expected Value of Total Tree Length")

ggplot(Data, aes(genecopies, TTL2)) +
  geom_point(color = "purple") +
  labs(x = "Number of Gene Copies", y = "Expected Number of Generations in N units", title = "Expected Value of Total Tree Length")
