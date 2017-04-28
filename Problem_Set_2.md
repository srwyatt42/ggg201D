# ProblemSet2
# Population Genetics Problem Set 2  
## Sydney Wyatt  
### May 1, 2017  


### Question 1

**A) Use a program such as R or Excel to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on initial allele frequency. The x-axis should be initial allele frequency and range from 0 to 1. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 0.1 for a population size of 2N=20.**  

```r
#Define function calculating expected delta f given initial allele f (fG1)
Expected.delta.f <- function(fG1) {
  #Create and empty data frame that you will eventually fill in
  delta.f.values = as.data.frame(cbind(rep(NA,21), rep(NA,21)))
  #Loop through all possible number of alleles given N = 10 and fill in the data frame
  for(i in 0:20){
    delta.f = i/20
    if(0< (fG1 + delta.f) & (fG1+delta.f) <=1){
      fG2 = fG1 + delta.f
      n = fG2*(2*N)
      P.deltaf.isX.a = dbinom(n, size = (2*N), prob = fG1)}
    else {
      P.deltaf.isX.a = 0}
    if(0<= (fG1 + delta.f) & (fG1+delta.f) <1){
      fG2 = fG1 - delta.f
      n = fG2*(2*N)
      P.deltaf.isX.b = dbinom(n, size = (2*N), prob = fG1)}
    else {
      P.deltaf.isX.b = 0}
    P.deltaf.isX = P.deltaf.isX.a + P.deltaf.isX.b
    delta.f.values[i+1,1] = delta.f
    delta.f.values[i+1,2] = P.deltaf.isX}
  #The data frame is full, now multiply the frequency (delta.f.values$V1) by the probability of X = x (delta.f.values$V2) and sum
  sum(delta.f.values$V1 * delta.f.values$V2)
}

N = 10
#Hold our values for the graph
Values.for.graph = as.data.frame(cbind(rep(NA, 11), rep(NA, 11)))

#Now, find the expected change in f for initial allele frequencies 0-1
for(i in 0:10){
  fG1 = i/10
  Values.for.graph[i+1, 1] = fG1
  Values.for.graph[i+1, 2] = Expected.delta.f(fG1)}

#Graph
plot(Values.for.graph, xlab = "Initial Allele Frequency", ylab = "Expected Change in Allele Frequency", main = "Expected Change in Allele Frequency Given Constant Population Size and Variable Initial Allele Frequency")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->




**B) Use the same program to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on population size. The x-axis should be population size and range from 2N=10 to 2N=100. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 10 with an allele frequency of 0.5.**  

```r
allele.freq = 0.5
initial.size = c(10,20,30,40,50,60,70,80,90,100)
pop.size.change = c(0,10,20,30,40,50,60,70,80,90,100)
Expected.change.pop = c()
```


### Question 2  

**A) Use the same program to generate a scatter plot that shows the properties of selection in large populations. The x-axis should be frequency of the advantageous allele and range from 0 to 1. The y-axis should be the change in frequency of the advantageous allele after one generation of selection. Perform calculations in steps of 0.01 for each of the following six (1a, 1b, 1c, 2a, 2b, 2c) scenarios: (1) the homozygous deleterious genotype has a selection coefficient of 0.1 and the advantageous allele is (a) recessive, (b) dominant or (c) additive; (2) the homozygous deleterious genotype has a selection coefficient of 0.25 and the advantageous allele is (a) recessive, (b) dominant or (c) additive.**  


```r
#Define function for recessive advantageous allele given S of "bad" genotype & frequency of "good" allele
recessive <- function(Sb, fg) {
  Sg = 0
  Shet = Sb
  fb = 1 - fg
  New.freq = ((fg^2)*(1-Sg)+(fg)*(fb)*(1-Shet))/((fg^2)*(1-Sg)+2*(fg)*(fb)*(1-Shet)+(fb^2)*(1-Sb))
  Change.in.freq = New.freq - fg
  return(Change.in.freq)
}

#Define function for dominant advantageous allele given S of "bad" genotype & frequency of "good" allele
dominant <- function(Sb, fg) {
  Sg = Shet = 0
  fb = 1 - fg
  New.freq = ((fg^2)*(1-Sg)+(fg)*(fb)*(1-Shet))/((fg^2)*(1-Sg)+2*(fg)*(fb)*(1-Shet)+(fb^2)*(1-Sb))
  Change.in.freq = New.freq - fg
  return(Change.in.freq)
}

#Define function for additive advantageous allele given S of "bad" genotype & frequency of "good" allele
#Assuming that selection for heterozygotes is half of the "bad" genotype's selection
additive <- function(Sb, fg) {
  Sg = 0
  Shet = 0.5*Sb
  fb = 1 - fg
  New.freq = ((fg^2)*(1-Sg)+(fg)*(fb)*(1-Shet))/((fg^2)*(1-Sg)+2*(fg)*(fb)*(1-Shet)+(fb^2)*(1-Sb))
  Change.in.freq = New.freq - fg
  return(Change.in.freq)
}
```


```r
#Now we use the above functions to calculate frequency of "good" allele in next generation

#Define allele frequency (This changes!)
good.allele.freq = seq(0,1, by = 0.01)   #0, 0.01, 0.02, ..... 1
bad.S = 0.1

#If "good" allele is recessive:
NextGenfreq.recessive = c()
for(i in good.allele.freq){
  f = recessive(bad.S, i)
  NextGenfreq.recessive = append(NextGenfreq.recessive, f)
}

#plot it!
plot(good.allele.freq, NextGenfreq.recessive, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "1a: S = 0.1, allele is recessive")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
#If "good" allele is dominant:
NextGenfreq.dominant = c()
for(i in good.allele.freq){
  f = dominant(bad.S, i)
  NextGenfreq.dominant = append(NextGenfreq.dominant, f)
}

plot(good.allele.freq, NextGenfreq.dominant, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "1b: S = 0.1, allele is dominant")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
#If "good" allele is additive:
NextGenfreq.additive = c()
for(i in good.allele.freq){
  f = additive(bad.S, i)
  NextGenfreq.additive = append(NextGenfreq.additive, f)
}

plot(good.allele.freq, NextGenfreq.additive, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "1c: S = 0.1, allele is additive")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

```r
#Now we have a new S!
bad.S = 0.25

#Copy the previous code and relabel the graphs.
#If "good" allele is recessive:
NextGenfreq.recessive = c()
for(i in good.allele.freq){   #iterate over the different "good" allele frequencies
  f = recessive(bad.S, i)
  NextGenfreq.recessive = append(NextGenfreq.recessive, f)
}

#plot it!
plot(good.allele.freq, NextGenfreq.recessive, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "2a: S = 0.25, allele is recessive")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

```r
#If "good" allele is dominant:
NextGenfreq.dominant = c()
for(i in good.allele.freq){
  f = dominant(bad.S, i)
  NextGenfreq.dominant = append(NextGenfreq.dominant, f)
}

plot(good.allele.freq, NextGenfreq.dominant, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "2b: S = 0.25, allele is dominant")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-5.png)<!-- -->

```r
#If "good" allele is additive:
NextGenfreq.additive = c()
for(i in good.allele.freq){
  f = additive(bad.S, i)
  NextGenfreq.additive = append(NextGenfreq.additive, f)
}

plot(good.allele.freq, NextGenfreq.additive, xlab = "Initial Allele Frequency", ylab = "Change in Allele Frequency", main = "2c: S = 0.25, allele is additive")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-6-6.png)<!-- -->





**B) Four of the six plots from above are highly asymmetric. Explain the biological reason behind these asymmetric patterns.**  



