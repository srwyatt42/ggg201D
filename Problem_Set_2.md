# ProblemSet2
# Population Genetics Problem Set 2  
## Sydney Wyatt  
### May 1, 2017  


### Question 1

**A) Use a program such as R or Excel to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on initial allele frequency. The x-axis should be initial allele frequency and range from 0 to 1. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 0.1 for a population size of 2N=20.**  

```r
initial.allele.freq = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
changing.allele.freq = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
N = 20
Expected.change = c()

for(i in initial.allele.freq){
  expected.change = c()
  for(p in changing.allele.freq){
    if (p==0){
      x = dbinom((20*i), size = N, prob = i)
      expected.change = append(expected.change, p*x)}
    else if (p==1){
      y = dbinom((20*p), size = N, prob = i)
      expected.change = append(expected.change,p*y)}
    else{
      z = dbinom((20*(i-0.1)), size = N, prob = i) + dbinom((20*(i+0.1)), size = N, prob = i)  
      expected.change = append(expected.change,p*z)}
  }
  Expected.change = append(Expected.change, (sum(expected.change)))
}

catplot(initial.allele.freq, Expected.change, size = 0.1, cat = 11, xlab = "Initial Allele Frequency", ylab = "New Allele Frequency", main = "Expected Change in Allele Frequency After One Generation")
```

![](Problem_Set_2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
## $xs
##  [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
## 
## $ys
##  [1] 0.0000000 0.9510997 1.1070411 1.1016797 1.0869896 1.0812101 1.0870262
##  [8] 1.1024776 1.1185703 1.0726763 1.0000000
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




**B) Use the same program to generate a scatter plot that shows how expected allele frequency change from genetic drift depends on population size. The x-axis should be population size and range from 2N=10 to 2N=100. The y-axis should be expected change in allele frequency after one generation. Perform calculations in steps of 10 with an allele frequency of 0.5.**  

```r
allele.freq = 0.5
initial.size = c(10,20,30,40,50,60,70,80,90,100)
pop.size.change = c(0,10,20,30,40,50,60,70,80,90,100)
Expected.change.pop = c()

for (i in initial.size){
  expected.change = c()
  for (n in pop.size.change){
    if (n==0){
      x = dbinom(i, size = i, prob = allele.freq)
      expected.change = append(expected.change, allele.freq*x)
    }
    if (n==100){
      y = dbinom(n, size = n, prob = allele.freq)
      expected.change = append(expected.change, allele.freq*y)
    }
    else {
      z = dbinom((i+n), size = i, prob = allele.freq) + dbinom((i-n), size = i, prob = allele.freq)
      expected.change = append(expected.change, allele.freq*z)
    }
  }
  Expected.change.pop = append(Expected.change.pop, sum(expected.change))
}
```


### Question 2  

**A) Use the same program to generate a scatter plot that shows the properties of selection in large populations. The x-axis should be frequency of the advantageous allele and range from 0 to 1. The y-axis should be the change in frequency of the advantageous allele after one generation of selection. Perform calculations in steps of 0.01 for each of the following six (1a, 1b, 1c, 2a, 2b, 2c) scenarios: (1) the homozygous deleterious genotype has a selection coefficient of 0.1 and the advantageous allele is (a) recessive, (b) dominant or (c) additive; (2) the homozygous deleterious genotype has a selection coefficient of 0.25 and the advantageous allele is (a) recessive, (b) dominant or (c) additive.**  



**B) Four of the six plots from above are highly asymmetric. Explain the biological reason behind these asymmetric patterns.**  


