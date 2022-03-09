#####splatter simulate
require(scater)
require(R.matlab)
source('./New_Metrics.R')
source('./DE_functions.R')
source('./Filter_functions.R')
source('./auroc.R')

# nts<-5
###simulated data selected

  dataset_name<-'EMTAB_h_k'
  
  data_new<-readRDS(paste0(dataset_name,'.rds'))
  datarow<-data_new[,which(colnames(data_new)==unique(colnames(data_new))[1])]
  #instances building & analysing
  ncells=dim(datarow)[2]
  maxsz<-floor(ncells*0.5)
  datarow<-datarow[which(rowSums(datarow)>0),]
  # dir.create(paste0('C:/R project/DEmethods/simulate buquan/',dataset_name,'/'))
  
  ######splatter

  
  for (j in c(0.2,0.4,0.6,0.8,1.0)) {
    
    #maxnum group 1 instance & analyse
    
    set.seed(100*j)
    idx_sample<-sample(ncol(datarow),2*floor(ncol(datarow)*j),replace = T)
    data_ori<-datarow[,idx_sample]
    group=c(rep('1',ncol(data_ori)/2),rep('2',ncol(data_ori)/2))
    group=factor(group)
    
    ### appeared later in group as cell type 2
    idx2<-which(group==unique(group)[2])
    ##cal_fc for each DE genes detected
    
    data_ori<-data_ori[which(rowSums(data_ori)>0),]
    data_ori<-as.matrix(data_ori)
    group_ori<-group
    
    DEmethods<-rev(c('scDD','MAST','t_test','BPSC','wilcox_test','samr','DESeq2','limma','DEsingle','edgeR'))
    filtermethods<-c(1:4)
   
    ###DE gene number detected by each methods record
    de_num_r=NULL
    
    cell_num_r<-NULL
    #######fpr
    fpr_record<-NULL
    ########TPR record
   
    for(filt in filtermethods){
      if(filt==1){
        ###filtering by OGFSC
        idx<-scCODE.filter_OGFSC(data_ori)
      }else{
        if(filt==2){
          ###filtering by conquer
          idx<-scCODE.filter_conquer(data_ori)
        }else{
          if(filt==3){
            ###filtering by scmap
            idx<-scCODE.filter_scmap(data_ori)
          }else{
            ###no filter
            idx=which(rowSums(data_ori)>0)
          }
        }
      }
      
      datatest<-data_ori[idx,]
      idx_c=which(colSums(datatest)==0)
      if(length(idx_c)!=0){
        datatest=datatest[,-idx_c]
        group=group_ori[-idx_c]
      }else{
        group=group_ori
      }
      datatest<-as.matrix(datatest)
      ###gene name
      genename<-rownames(datatest)
      
      for (de in DEmethods) {
        
        
        run_func<-paste0('scCODE.',de,'(datatest,group)')
        res_temp<-try(eval(parse(text = run_func)),silent = T)
        if('try-error'%in%class(res_temp)){
          res_temp=NA
        }
        ### adjust p value<=0.05
       gene_num<-nrow(datatest)
       de_num<-length(which(res_temp<=.05))
       fpr_temp<-de_num/gene_num
       
       cell_num_r<-c(cell_num_r,ncol(datatest)/2)
       de_num_r<-c(de_num_r,de_num)
       
       fpr_record<-c(fpr_record,fpr_temp)
        
      }
      
    }
   
    dir.create(paste0('./',dataset_name))
    fpr_null_on_simulate<-cbind(cell_num_r,de_num_r,fpr_record)
    rownames(fpr_null_on_simulate)<-c(paste0('OGFSC-',DEmethods),paste0('conquer-',DEmethods),paste0('scmap-',DEmethods),paste0('without_filter-',DEmethods))
    write.csv(fpr_null_on_simulate,paste0('./',dataset_name,'/Evaluate-fpr_null-',j,'.csv'))
    
    cat(paste0('/nAnalysis Finished/n dep-'),j)
  }
  

