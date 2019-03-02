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

### 2019.2.22
kBET exploration

### 2019.3.2
quiz-ROC plot

quiz.ipynb (use diff_exp feature to plot ROC, fpr-tpr dataframe provided)
quiz_tired.ipynb to make feature-selection files used for ROC plot
quiz_right.ipynb (ROC plot

