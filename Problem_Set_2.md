# ProblemSet2
# Population Genetics Problem Set 2  
## Sydney Wyatt  
### May 1, 2017  


### Question 1

**A) Use a program such as R or Excel to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on initial allele frequency. The x-axis should be initial allele frequency and range from 0 to 1. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 0.1 for a population size of 2N=20.**  

```r
A.initialprob = seq(0.1,1, by = 0.1)
A.changingprob = seq(0,1, by = 0.1)
A.n = 20
A.SecondGenExpChange = c()

for(i in A.initialprob){
  expected.change = c()
  for(p in seq_along(A.changingprob)){
    if (p==1){
      x = A.changingprob[p]*(dbinom((20*A.changingprob[p]), size = A.n, prob = i) + dbinom((20*A.changingprob[p]+1), size = A.n, prob = i))
      expected.change = append(expected.change, x)}
    else if (p==11){
      y = A.changingprob[p]*(dbinom((20*A.changingprob[p]-1), size = A.n, prob = i) + dbinom((20*A.changingprob[p]), size = A.n, prob = i))
      expected.change = append(expected.change,y)}
    else{
      z = A.changingprob[p]*(dbinom((20*A.changingprob[p]-1), size = A.n, prob = i) + dbinom((20*A.changingprob[p]+1), size = A.n, prob = i))  
      expected.change = append(expected.change,z)}
  }
  A.SecondGenExpChange = append(A.SecondGenExpChange, sum(expected.change))
}

catplot(A.initialprob, A.SecondGenExpChange, size = 0.1, cat = 11, xlab = "Initial Allele Frequency", ylab = "New Allele Frequency", main = "Expected Change in Allele Frequency After One Generation")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
## $xs
##  [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
## 
## $ys
##  [1] 0.1014412 0.2000122 0.3000000 0.4000000 0.5000010 0.6000366 0.7007979
##  [8] 0.8114805 1.0086063 1.0000000
## 
## $args
## $args$xlab
## [1] "Initial Allele Frequency"
## 
## $args$ylab
## [1] "New Allele Frequency"
## 
## $args$main
## [1] "Expected Change in Allele Frequency After One Generation"
## 
## 
## $canvas
## [1] 0.0 1.1 0.0 1.1
```

```r
#Simple for and if else loops:
#Not.Div = c()
#Div = c()
#for(i in 1:10){
#  if (i %% 4){
#    Not.Div = append(Not.Div,i)}
#  else {
#    Div = append(Div, i)}
#}
```


**B) Use the same program to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on population size. The x-axis should be population size and range from 2N=10 to 2N=100. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 10 with an allele frequency of 0.5.**  



### Question 2  

**A) Use the same program to generate a scatter plot that shows the properties of selection in large populations. The x-axis should be frequency of the advantageous allele and range from 0 to 1. The y-axis should be the change in frequency of the advantageous allele after one generation of selection. Perform calculations in steps of 0.01 for each of the following six (1a, 1b, 1c, 2a, 2b, 2c) scenarios: (1) the homozygous deleterious genotype has a selection coefficient of 0.1 and the advantageous allele is (a) recessive, (b) dominant or (c) additive; (2) the homozygous deleterious genotype has a selection coefficient of 0.25 and the advantageous allele is (a) recessive, (b) dominant or (c) additive.**  



**B) Four of the six plots from above are highly asymmetric. Explain the biological reason behind these asymmetric patterns.**  

