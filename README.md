# exSeek_package
There are some tiny probelms in the scripts Zhuoer wrote for matrix processing. Imputation and batch removal were not included in matrix processing. 
## WHAT I DO IS:
* Add imputation.
* Add batch removal including RUVs and combat.
* Remove logic error.
* ...

[rexSeek](https://github.com/dongzhuoer/rexseek)

### 2018.12.03
|package|author|
| :-- | :-- |
|normalization.R|zhuoer|
|matrix-process.R|zhuoer|
|batch-removal.R|yufengxie|

#### I improved scripts zhuoer wrote,
* matrix-process.R

add imputation function

* normalization.R
main function zhuoer wrote lacking imputation, I add this part into normalization main function.

#### I wrote script for batch removal,
* batch-removal.R
##### reference:
```
origin script:
home/chenxupeng/projects/exseek/jupyter/matrix_processing.ipynb

combat:
http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html
https://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
```

