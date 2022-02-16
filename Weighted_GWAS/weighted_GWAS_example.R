#Weighted GWAS using example genotypes at 1000 SNPs - n.b. sample order has been shuffled
#Code to adjust standard error for reduced effective sample size (lines: 22-23) courtesy of James Watson

#read in sample genotypes
geno <- read.table("example_geno_chunk.txt.gz", header = F)
#read in covariate, case/control, weights data
samples <- read.table("samples.txt", header = T)


geno.t <- t(geno[,c(7:dim(geno)[2])])
total <- cbind(samples, geno.t)

#loop weighted logistic regression association analysis over SNPs
pval <- c()
se <- c()
beta <- c()

for (i in c(10:dim(total)[2])){
  
  fit_glm1 = glm(total$cc ~ total[,i] + total$PC1 + total$PC2 + total$PC3 + total$PC4 + total$PC5 + total$PC6 + total$PLATFORM, family='quasibinomial', 
                 weights=total$weights)
  
  cov.m1 <- sandwich::vcovHC(fit_glm1, type = "HC0")
  std.err <- sqrt(diag(cov.m1))[2]
  beta[(i-9)] <- coef(fit_glm1)[2]
  pval[(i-9)] <- 2 * pnorm(abs(coef(fit_glm1)[2]/std.err), lower.tail = FALSE)
  se[(i-9)] <- std.err
}

#results
res.out <- cbind(geno[,c(1,4,3,5,6)], pval, beta, se)

