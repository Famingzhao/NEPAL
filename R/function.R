#' Title NEPAL_bulk
#'
#' @param bulk.data input bulk transcriptome data.
#' @param method
#' @param species
#' @param gene.sets
#'
#' @return
#' @export
#'
#' @examples
NEPAL_bulk <- function(bulk.data,
                             method=c("ssGSEA","Enet","Ridge","SVM","all"),
                             species=c("human","mouse"),
                             gene.sets = c("all","EurUrol.2005", "CancerDiscov.2011", "CancerRes.2014", "NatMed.2016",
                                           "CellRep.2018.scRNA", "ClinCancerRes.2015", "JClinInvest.2019",
                                           "IntJCancer.2019", "HP_NE_neoplasm", "JClinOncol.2018", "BMC.Cancer.2017",
                                           "NEPAL")
){
  species <- match.arg(species)
  gene.sets <- match.arg(gene.sets)
  suppressWarnings(library(dplyr))
  suppressWarnings(library(glmnet))
  suppressWarnings(library(e1071))
  suppressWarnings(library(IOBR))
  seed <- 1
  load(system.file("data", "model_gene.sets.Rdata",
                   package = "NEPAL", mustWork = TRUE))
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }

  if(species == "mouse"){
    load(system.file("data", "biomaRt_ensembl_annotation.Rdata",
                     package = "NEPAL", mustWork = TRUE))
    convertMouseGeneList.Bulk <- function(bulk.data,
                                          Mus.data = F,
                                          Human.data = F,...){
      count.data = bulk.data

      # mouse转Human
      convertMouseGeneList <- function(x = x){
        require("biomaRt")
        # human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # save(human,mouse,file = "./biomaRt_ensembl_annotation.Rdata")
        genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                         values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        # humanx <- unique(genesV2[, 2])
        # Print the first 6 genes found to the screen
        # print(head(genesV2))
        return(genesV2)
      }

      # Human转mouse
      convertHumanGeneList <- function(x){
        require("biomaRt")
        # human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # save(human,mouse,file = "./biomaRt_ensembl_annotation.Rdata")
        genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                         values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
        # humanx <- genesV2
        # Print the first 6 genes found to the screen
        # print(head(genesV2))
        return(genesV2)
      }

      if(Mus.data){
        hum_genes <- convertMouseGeneList(x = row.names(count.data))
        head(hum_genes)
        hum_genes = hum_genes[!is.na(hum_genes$HGNC.symbol)|hum_genes$HGNC.symbol=="",]

        #### merge
        count.data = as.data.frame(count.data) %>% rownames_to_column(var = "MGI.symbol")
        count.data = left_join(count.data,hum_genes)

        count.data = count.data %>%
          dplyr::select(HGNC.symbol,MGI.symbol,everything())
        count.data = count.data[!is.na(count.data$HGNC.symbol),]

        rank.data = rowSums(count.data[,-c(1,2)])

        count.data = count.data[order(rank.data,decreasing = T),]
        count.data = count.data[!duplicated(count.data$HGNC.symbol),]
        row.names(count.data) =  count.data$HGNC.symbol
        count.data = count.data[,-c(1,2)]
      }

      if(Human.data){
        mus_genes <- convertHumanGeneList(x = row.names(count.data))
        head(mus_genes)
        mus_genes = mus_genes[!is.na(mus_genes$MGI.symbol)|mus_genes$MGI.symbol=="",]

        #### 3.merge
        count.data = as.data.frame(count.data) %>% rownames_to_column(var = "HGNC.symbol")
        count.data = left_join(count.data,mus_genes)

        count.data = count.data %>%
          dplyr::select(HGNC.symbol,MGI.symbol,everything())

        count.data = count.data[!is.na(count.data$MGI.symbol),]

        rank.data = rowSums(count.data[,-c(1,2)])

        count.data = count.data[order(rank.data,decreasing = T),]
        count.data = count.data[!duplicated(count.data$MGI.symbol),]
        row.names(count.data) =  count.data$MGI.symbol
        count.data = count.data[,-c(1,2)]
      }
      return(count.data)
    }
    bulk.data = convertMouseGeneList.Bulk(bulk.data = bulk.data,Mus.data = T,Human.data = F)
    input.data = bulk.data[cand.genes,]
  }
  # input data
  if(species == "human"){
    input.data = bulk.data[cand.genes,]
  }
  row.names(input.data) = cand.genes
  input.data[is.na(input.data)] = 0

  input.data <- t(input.data) %>% as.data.frame()
  input.data <- input.data %>% tibble::rownames_to_column("ID")
  val_dd_list = list(input_data = input.data)
  ##################################
  #### Elastic net ####
  ##################################
  if(method=="Enet"){
    rs = lapply(names(Enet.res.list), function(alpha){
      set.seed(seed)
      rs <- lapply(val_dd_list, function(x){
        cbind(x[,1,drop=F],
              RS = as.numeric(predict(Enet.res.list[[alpha]],type='link',
                                      newx=as.matrix(x[,-c(1)]),
                                      s=Enet.res.list[[alpha]]$lambda.min)) %>% normalize())})
      rs = lapply(rs, function(x){
        colnames(x)[1:2] = c("ID",alpha)
        return(x)
      })
      return(rs)
    })
    names(rs) = names(Enet.res.list)

    finnal.res.list = lapply(names(val_dd_list), function(data.sets){
      fin.res = val_dd_list[[data.sets]][,1,drop=F]
      for (alpha in names(rs)) {
        test = rs[[alpha]][[data.sets]]
        fin.res =  left_join(fin.res, test)
      }
      return(fin.res)
    })
    names(finnal.res.list) = names(val_dd_list)
    rs = finnal.res.list$input_data
  }

  ##################################
  #### Ridge ####
  ##################################
  if(method=="Ridge"){
    rs <- lapply(val_dd_list, function(x){
      cbind(x[,1,drop=F],
            Ridge.score = as.numeric(predict(Ridge.fit,type='response',
                                             newx=as.matrix(x[,-c(1)]),
                                             s=Ridge.fit$lambda.min)) %>% normalize())})
    rs = rs$input_data
  }

  ##################################
  #### SVM ####
  ##################################
  if(method=="SVM"){
    rs <- lapply(val_dd_list, function(x){
      cbind(x[,1,drop=F],
            SVM.score = as.numeric(predict(fit.svm,x)) %>% normalize())})
    rs = rs$input_data
  }

  ##################################
  #### ssGSEA ####
  #################################
  if(method=="ssGSEA"){
    if(gene.sets=="all"){
      rs <- IOBR::calculate_sig_score(pdata           = NULL,
                                      eset            = bulk.data,
                                      signature       = ne.markers.list,
                                      method          = "ssgsea",
                                      mini_gene_count = 0)
      rs = arrange(rs,Index)
      rs$NE_UP_DN = rs$NE_sig_UP-rs$NE_sig_DN
    }

    if(gene.sets=="NEPAL"){
      rs <- IOBR::calculate_sig_score(pdata           = NULL,
                                      eset            = bulk.data,
                                      signature       = ne.markers.list[c("NE_sig_UP",
                                                                          "NE_sig_DN")],
                                      method          = "ssgsea",
                                      mini_gene_count = 0)
      rs = arrange(rs,Index)
      rs$NE_UP_DN = rs$NE_sig_UP-rs$NE_sig_DN
    } else{
      rs <- IOBR::calculate_sig_score(pdata           = NULL,
                                      eset            = bulk.data,
                                      signature       = ne.markers.list[gene.sets],
                                      method          = "ssgsea",
                                      mini_gene_count = 0)
      rs = arrange(rs,Index)
    }
  }

  ##################################
  #### All ####
  ##################################
  if(method == "all"){
    #Enet
    rs.Enet.list = lapply(names(Enet.res.list), function(alpha){
      set.seed(seed)
      rs <- lapply(val_dd_list, function(x){
        cbind(x[,1,drop=F],
              RS = as.numeric(predict(Enet.res.list[[alpha]],type='link',
                                      newx=as.matrix(x[,-c(1)]),
                                      s=Enet.res.list[[alpha]]$lambda.min)) %>% normalize())})
      rs = lapply(rs, function(x){
        colnames(x)[1:2] = c("ID",alpha)
        return(x)
      })
      return(rs)
    })
    names(rs.Enet.list) = names(Enet.res.list)

    finnal.res.list = lapply(names(val_dd_list), function(data.sets){
      fin.res = val_dd_list[[data.sets]][,1,drop=F]
      for (alpha in names(rs.Enet.list)) {
        test = rs.Enet.list[[alpha]][[data.sets]]
        fin.res =  left_join(fin.res, test)
      }
      return(fin.res)
    })
    names(finnal.res.list) = names(val_dd_list)
    rs.Enet.list = finnal.res.list
    str(rs.Enet.list)

    #Ridge
    rs.Ridge <- lapply(val_dd_list, function(x){
      cbind(x[,1,drop=F],
            Ridge.score = as.numeric(predict(Ridge.fit,type='response',
                                             newx=as.matrix(x[,-c(1)]),
                                             s=Ridge.fit$lambda.min)) %>% normalize())})

    #SVM
    rs.SVM <- lapply(val_dd_list, function(x){
      cbind(x[,1,drop=F],
            SVM.score = as.numeric(predict(fit.svm,x)) %>% normalize())})

    #ssGSEA
    res.ssGSEA <- IOBR::calculate_sig_score(pdata           = NULL,
                                            eset            = bulk.data,
                                            signature       = ne.markers.list[c("NE_sig_UP",
                                                                                "NE_sig_DN")],
                                            method          = "ssgsea",
                                            mini_gene_count = 0)
    res.ssGSEA = arrange(res.ssGSEA,Index)
    res.ssGSEA$NE_UP_DN = res.ssGSEA$NE_sig_UP-res.ssGSEA$NE_sig_DN

    # merge
    res.sum = lapply(names(val_dd_list), function(data.sets){
      data = left_join(res.ssGSEA,rs.Enet.list[[data.sets]]) %>%
        left_join(rs.Ridge[[data.sets]]) %>% left_join(rs.SVM[[data.sets]])
      return(data)
    })
    rs = res.sum[[1]]
  }
  return(as.data.frame(rs))
}


#' Title NEPAL_scRNA
#'
#' @param seurat.data input scRNA-seq seurat data.
#' @param method
#' @param species
#' @param ncores
#' @param assay.names
#' @param DefaultAssay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
NEPAL_scRNA <- function(seurat.data,
                              method = "AUCell",
                              species=c("human","mouse"),
                              ncores=2,assay.names = "NE",DefaultAssay=F,...
){
  species <- match.arg(species)
  suppressWarnings(library(dplyr))
  sc.Pathway.Seurat = function (obj, method = "AUCell", imputation = F, ncores = 2,
                                geneList = system.file("data", "NE_cells_markers.gmt",package = "NEPAL"),
                                assay.names = "NE",DefaultAssay=F)

  {
    gmtFile = geneList
    geneSets <- GSEABase::getGmt(geneList)
    ##########
    countexp <- obj@assays$RNA@counts
    countexp <- data.frame(as.matrix(countexp))
    ##########
    if (imputation == F) {
      countexp2 <- countexp
    }
    if (imputation == T) {
      library(rsvd)
      cat("Start imputation...\n")
      cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
      result.completed <- alra(as.matrix(t(countexp)))
      countexp2 <- result.completed[[3]]
      row.names(countexp2) <- colnames(countexp)
    }
    cat("Start quantify the pathway activity...\n")
    if (method == "VISION") {
      library(VISION)
      n.umi <- colSums(countexp2)
      scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
      vis <- Vision(scaled_counts, signatures = gmtFile)
      options(mc.cores = ncores)
      vis <- analyze(vis)
      signature_exp <- data.frame(t(vis@SigScores))
    }
    if (method == "AUCell") {
      library(AUCell)
      library(GSEABase)
      cells_rankings <- AUCell_buildRankings(as.matrix(countexp2),
                                             nCores = ncores, plotStats = F)
      geneSets <- getGmt(gmtFile)
      cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
      signature_exp <- data.frame(getAUC(cells_AUC))
      signature_exp["NE_UP_DN",] = signature_exp["NE_sig_UP",] - signature_exp["NE_sig_DN",]
      signature_exp["NE_UP_DN",] = ifelse(as.numeric(signature_exp["NE_UP_DN",])<0,0,
                                          as.numeric(signature_exp["NE_UP_DN",]))
    }
    if (method == "ssGSEA") {
      library(GSVA)
      library(GSEABase)
      geneSets <- getGmt(gmtFile)
      gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"),
                      kcdf = c("Poisson"), parallel.sz = ncores)
      signature_exp <- data.frame(gsva_es)
    }
    if (method == "gsva") {
      library(GSVA)
      library(GSEABase)
      geneSets <- getGmt(gmtFile)
      gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"),
                      kcdf = c("Poisson"), parallel.sz = ncores)
      signature_exp <- data.frame(gsva_es)
    }
    colnames(signature_exp) = row.names(obj@meta.data)
    obj[[assay.names]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
    obj <- SeuratObject::SetAssayData(obj, slot = "scale.data",
                                      new.data = as.matrix(signature_exp), assay = assay.names)
    if(DefaultAssay == T){
      DefaultAssay(obj) <- assay.names
    }
    obj
  }
  if(species == "mouse"){
    load(system.file("data", "biomaRt_ensembl_annotation.Rdata",
                     package = "NEPAL", mustWork = TRUE))
    convertMouseGeneList.Seurat <- function(seurat.data,
                                            Mus.data = F,
                                            Human.data = F,...){
      count.data = seurat.data@assays$RNA@counts
      meta.data = seurat.data@meta.data

      # mouse转Human
      convertMouseGeneList <- function(x = x){
        require("biomaRt")
        # human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # save(human,mouse,file = "./biomaRt_ensembl_annotation.Rdata")
        genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                         values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        # humanx <- unique(genesV2[, 2])
        # Print the first 6 genes found to the screen
        # print(head(genesV2))
        return(genesV2)
      }

      # Human转mouse
      convertHumanGeneList <- function(x){
        require("biomaRt")
        # human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
        # save(human,mouse,file = "./biomaRt_ensembl_annotation.Rdata")
        genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                         values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
        # humanx <- genesV2
        # Print the first 6 genes found to the screen
        # print(head(genesV2))
        return(genesV2)
      }

      if(Mus.data){
        hum_genes <- convertMouseGeneList(x = row.names(count.data))
        head(hum_genes)
        hum_genes = hum_genes[!is.na(hum_genes$HGNC.symbol)|hum_genes$HGNC.symbol=="",]

        #### merge
        count.data = as.data.frame(count.data) %>% rownames_to_column(var = "MGI.symbol")
        count.data = left_join(count.data,hum_genes)

        count.data = count.data %>%
          dplyr::select(HGNC.symbol,MGI.symbol,everything())
        count.data = count.data[!is.na(count.data$HGNC.symbol),]

        rank.data = rowSums(count.data[,-c(1,2)])

        count.data = count.data[order(rank.data,decreasing = T),]
        count.data = count.data[!duplicated(count.data$HGNC.symbol),]
        row.names(count.data) =  count.data$HGNC.symbol
        count.data = count.data[,-c(1,2)]
        sec = CreateSeuratObject(counts = count.data,
                                 meta.data = meta.data,
                                 min.cells = 0,
                                 min.features = 0)
      }

      if(Human.data){
        mus_genes <- convertHumanGeneList(x = row.names(count.data))
        head(mus_genes)
        mus_genes = mus_genes[!is.na(mus_genes$MGI.symbol)|mus_genes$MGI.symbol=="",]

        #### 3.merge
        count.data = as.data.frame(count.data) %>% rownames_to_column(var = "HGNC.symbol")
        count.data = left_join(count.data,mus_genes)

        count.data = count.data %>%
          dplyr::select(HGNC.symbol,MGI.symbol,everything())

        count.data = count.data[!is.na(count.data$MGI.symbol),]

        rank.data = rowSums(count.data[,-c(1,2)])

        count.data = count.data[order(rank.data,decreasing = T),]
        count.data = count.data[!duplicated(count.data$MGI.symbol),]
        row.names(count.data) =  count.data$MGI.symbol
        count.data = count.data[,-c(1,2)]
        sec = CreateSeuratObject(counts = count.data,
                                 meta.data = meta.data,
                                 min.cells = 0,
                                 min.features = 0)
      }
      return(sec)
    }
    seurat.data = convertMouseGeneList.Seurat(seurat.data = seurat.data,
                                              Mus.data = T,Human.data = F)
  }
  seurat.data = sc.Pathway.Seurat(obj = seurat.data,
                                  method = method,
                                  ncores = ncores,
                                  assay.names = assay.names,
                                  DefaultAssay=DefaultAssay)

  return(seurat.data)
}


#' Title draw.survival
#'
#' @param survfit
#' @param ylab
#' @param legend.position
#' @param plot.title
#' @param break.time
#' @param legend.labs
#' @param xlimt
#' @param risk.table
#' @param front.size
#' @param palette
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
draw.survival <- function(survfit= fit_PFI,
                          ylab= "BCR",
                          legend.position = "right",
                          plot.title = NULL,
                          break.time=NA,
                          legend.labs = c('NE.sig High','NE.sig Low'),
                          xlimt = NA,
                          risk.table = F,
                          front.size = 10,
                          palette = c("#1F78B4","#E31A1C"),...) {
  if(ylab %in% c("OS","RFS")){
    if(ylab== "OS"){
      ylab <- "Overall survival (%)"}
    if(ylab== "BCR"){
      ylab <- "BCR-free survival (%)"}
  }else{ylab = ylab}

  if(is.na(xlimt)){
    xlimt = ceiling(max(na.omit(survfit$time)))
  }else{ylab = ylab}

  if(is.na(break.time)){
    if(xlimt<80|xlimt==80){
      break.time = 20
    }else if(xlimt>80&xlimt<120 | xlimt==120){
      break.time = 30
    }else if(xlimt>120&xlimt<300){
      break.time = 50
    }else{break.time = 100}
  }else{ylab = ylab}

  if(T){mytheme <- theme(plot.title = element_text(size = front.size+2,color="black",hjust = 0.5),
                         axis.title = element_text(size = front.size,color ="black"),
                         axis.text = element_text(size=front.size,color = "black"),
                         #axis.line = element_line(color = "black"),
                         #axis.ticks = element_line(color = "black"),
                         #panel.grid.minor.y = element_blank(),
                         #panel.grid.minor.x = element_blank(),
                         # panel.grid=element_blank(),
                         legend.position = legend.position,
                         legend.text = element_text(size=front.size),
                         legend.title= element_text(size= front.size)
  )}
  pp_OS <- ggsurvplot(
    survfit,
    pval = TRUE,
    pval.method = TRUE,
    legend.title="",
    legend.labs= legend.labs,
    risk.table = risk.table,
    risk.table.y.text.col = risk.table,
    risk.table.y.text = risk.table,
    ggtheme = theme_bw(),
    palette = palette,
    xlab = 'Time in months',
    ylab=ylab,
    break.time.by = break.time )

  pp_OS$plot <- pp_OS$plot + labs(title = plot.title)+
    scale_y_continuous(labels = scales::percent)+
    coord_cartesian(xlim = c(0,xlimt)) + mytheme
  # 改变Risk Table 的label可以用pp$label
  if(risk.table){
    pp_OS$table <- pp_OS$table+labs(title = "Under risk",x="Time in months")+
      # theme(plot.title = element_text(size = 11,color="black"))+
      coord_cartesian(xlim = c(0,xlimt))+ mytheme
  }


  return(pp_OS)
}

