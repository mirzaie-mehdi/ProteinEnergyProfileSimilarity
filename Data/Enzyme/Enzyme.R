# Enzyme
df1 <- read.delim('Data/Enzyme/Hydrolases_Positive_Train.Txt')
df2 <- read.delim('Data/Enzyme/Isomerases_Positive_Train.txt')
df3 <- read.delim('Data/Enzyme/Ligases_Positive_Train.txt')
df4 <- read.delim('Data/Enzyme/Lyases_Positive_Train.txt')
df5 <- read.delim('Data/Enzyme/Oxidoreductases_Positive_Train.txt')
df6 <- read.delim('Data/Enzyme/Transferases_Positive_Train.txt')

size.class <- c(nrow(df1),nrow(df2),nrow(df3),nrow(df4),nrow(df5),nrow(df6))

E210.seq_Enzyme <- data.frame(class=rep(c('Hydrolases','Isomerases',
                                          'Ligases','Lyases',
                                          'Oxidoreductases','Transferases'),
                                        each=1000),
                              uniID=c(df1[sample(1:nrow(df1),1000),1],  #
                                      df2[sample(1:nrow(df2),1000),1],  #
                                      df3[sample(1:nrow(df3),1000),1],  #
                                      df4[sample(1:nrow(df4),1000),1],  #
                                      df5[sample(1:nrow(df5),1000),1],  #
                                      df6[sample(1:nrow(df6),1000),1]), #
                              length=NA,
                              matrix(0,6000,210))

Nprotein <- nrow(E210.seq_Enzyme)
#----------------------------------------------
colnames(E210.seq_Enzyme)[-c(1:3)] <- pair_inter

for (NP in 1:Nprotein) {
  print(NP)
  tryCatch({
    seq <- getUniProt(E210.seq_Enzyme$uniID[NP])[[1]]
    seq <- unlist(strsplit(seq,split = ''))
    E210.seq_Enzyme$length[NP] <- length(seq)
    freq <- c()
    for (j in 1:20){
      freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
    }
    freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
    en_fr <- aaenergy * freq_matrix
    en_fr <- en_fr/length(seq)
    pair_es <- diag(freq) %*% as.matrix(en_fr)
    aa210_p <- pair_es[lower.tri(pair_es,diag = T)]
    E210.seq_Enzyme[NP,4:213] <- aa210_p
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}
saveRDS(E210.seq_Enzyme,'Data/Enzyme/E210.seq_Enzyme.rds')
# -----------------------------------------
# -----------------  UMAP
source('mytheme.R')
E210.seq_Enzyme <- readRDS('Data/Enzyme/E210.seq_Enzyme.rds')
df <- E210.seq_Enzyme
df <- df[!is.na(df$length),]
dis <- as.matrix(proxy::dist(df2[,-c(1:3)], method = manhat))
ump <- umap(d = dis,
            n_neighbors = 100,
            min_dist = 0.1,
            input="dist")

ump_df <- data.frame(class=df2$class,
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

ggplot(ump_df,
           aes(UMAP1, UMAP2, color= class)) +
  geom_point(alpha=.7, size=2) + 
  mytheme()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom',
        axis.title = element_text(size = 15),
        axis.title.y = element_text(size = 15,vjust = 5),
        legend.text = element_text(size = 15))

# ---------------------------------------
# ----------- RandomForest
df$class <- factor(df$class)
kfolds <- createFolds(df$class, k = 10)
res<-c()
for (i in 1:10) {
  print(i)
  x <- kfolds[[i]]
  tr <- df[-x,-c(2:3)]
  te <- df[x,-c(2:3)]
  rf = randomForest(class~.,data = tr)
  p2 <- predict(rf, te)
  cm2 <- confusionMatrix(p2, te$superfamily)
  Accuracy <- c(cm2$overall['Accuracy'], F1=cm2$byClass[,'F1'])
  res<-rbind(res,Accuracy)
}
res <- rbind(res,colMeans(res))
res <- round(res,2)
rownames(res) <- c(paste0("Fold_",1:length(kfolds)),'Average')
colnames(res)[-1] <- gsub('.Class','',colnames(res)[-1])
# -----------------
tab1 <- rbind(res1nn,res[11,])
colnames(tab1)<-colnames(res)
tab1 <- data.frame(method=c('1NN','RF'),tab1)
