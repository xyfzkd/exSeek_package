# exSeek_package
There are some tiny probelms in the scripts Zhuoer wrote for matrix processing. Imputation and batch removal were not included in matrix processing.
## WHAT I DO IS:
* Add imputation.
* Add batch removal including RUVs and combat.
* Remove logic error.
* ...

[rexSeek](https://github.com/dongzhuoer/rexseek)

### 2018.12.03
|script|author|
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
### 2018.12.06
* plot.R

run well on my account, but failed in Binbin & Xupeng's

### 2018.12.08
* refer_output.ipynb

for function `norm_cpm_refer`, I add a command, allowing producing refer_id.txt

### 2018.12.10
top k function changed

top k genes use counts top k sum as factor
others use count down sum as factor

But what happened to counts value near cut-off? How to solve this problem?

### 2018.12.27
wrote plot codes of {R,python}

based of `plot.ipynb`, I wrote R version `r2py_plot-Copy1.ipynb` and python version `r2py_plot.ipynb`, the latter of which output bad pics because of rpy2 package.

#### some problems solved
* plot_highest_exprs using gene ID as index

***
diff_exp added
function diff_exp in `diff_exp.ipynb` section 6.2 can produce files such as `Healthy.CRC.csv` containing padj_values
[Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory)
[Package ‘DESeq2’](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)
