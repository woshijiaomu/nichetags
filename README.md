# nichetags install guide
you can install it from github:

```
install.packages(c("Polychrome","stringdist","stringr","ggplot2","cowplot","igraph"))
if(!require(devtools)){
	install.packages("devtools")
}
if(!require(nichetags)){
	library(devtools)
	install_github("woshijiaomu/nichetags")
}
```

after installation，attach the package in R：

```
library(nichetags)
```
