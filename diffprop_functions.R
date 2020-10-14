# Differential proportion analysis of single-cell populations

###
# Simulation of a data set with different number of cells and their cell types
###
simCells<-function(n=10^7, prop=c(0.6, 0.3, 0.05, 0.03, 0.01, 0.005, 0.005)){
	ct<-paste("Cell_Type_",1:length(prop),sep="");
	nct<-round(n*prop);
	lab<-rep(ct,nct);	
	return(lab);
}

###
# Simulate a data-set with some provided error rate to simulate experimental/biological noise
###
simCellsNoisy<-function(n=10^7, prop, error = 0.05){
  # Create a 'noisy' distribution for sampling cells
  prop.noisy = unlist(lapply(prop, function(x) 
                     {ifelse(sample(c(1:2), size=1) == 1, x+x*error, x-x*error)}))
  prop.noisy = prop.noisy/sum(prop.noisy)
  
  # Sample as before but with a modified distribution
  ct<-paste("Cell_Type_",1:length(prop.noisy),sep="");
  nct<-round(n*prop.noisy);
  lab<-rep(ct,nct);	
  return(lab);
}

###
# Subsampling of cells
###
subCells<-function(cells, n=round(length(cells)/10)){
	sample(cells,n);
}



###
# Create the condition and cluster variables from the observed contingency table
###
makeVariables<-function(obs.table){
   cond_names<-dimnames(obs.table)[[1]]; #index names
   clus_names<-dimnames(obs.table)[[2]]; #col names
   
   cond<-rep(cond_names, apply(obs.table,1,sum) ); #generates a list with each condition (treatment) lenght is totalnumber of cells
   clus_list<-apply(obs.table,1,function(x){ #list with cluster for each treatment 
     rep(clus_names, x ); # repeats something x times, there is some weird bug that sometimes the clus_list becomes a matrix, this is dependent on how it was normalized - don't fully understand what is going on...
   })
   clus<-base::do.call(c,clus_list); # make a table with condition (named ct1, ct2 in one column, cluster ID in other column)
   #table(cond,clus);
   return(list(cond=cond, clus=clus));
}

###
# Generate null distributions based on a permutation procedure
###

# n = number of permutations
# p = proportion of cells that are being permuted

generateNull<-function(obs.table, n=10000, p=0.2){
    obs<-makeVariables(obs.table);
    all.exp<-array(NA,dim=c(dim(obs.table)[1],dim(obs.table)[2],n)) #make a list that is 10000 long
    dimnames(all.exp)[[1]]<-dimnames(obs.table)[[1]];
    dimnames(all.exp)[[2]]<-dimnames(obs.table)[[2]];
    dimnames(all.exp)[[3]]<-1:n; # set the dimensions of the all.exp with cluster, treatment and 10000 replicstes
    clus_names<-dimnames(obs.table)[[2]];
    cond_names<-dimnames(obs.table)[[1]];
    
    v<-makeVariables(obs.table); # make the variables for each cell which treatment does it belong to and which clusters
    
    ### Perform permutation n times
    for(i in 1:n){
      pv<-v;
      
      if (i %% 1000 == 0){
        print(i);}
      
      ### Permutation: Randomly select p*100% of points to be re-sampled from the background 
        
      # this selects the number of cells to resample - this is based on the total number of conditions so 20% (=p) of total # of cells
      # randInd is a vector of cell indices (lenght is 20% of total number of cells, basically the new indices of cells)
      randInd<-sample(1:length(pv$clus),round(length(pv$clus)*p)); 
      
      # replace the indices of the sample lables - that's how sample labels are switched
      pv$clus[randInd]<-sample(v$clus,length(randInd),replace=F);  
      
      this.exp<-table(pv$cond,pv$clus);
      exp.table<-obs.table;
      exp.table[1:dim(exp.table)[1],1:dim(exp.table)[2]]<-NA
      exp.table[dimnames(this.exp)[[1]],dimnames(this.exp)[[2]]]<-this.exp;
      exp.table
      all.exp[,,i]<-exp.table;
    }
    return(all.exp);
}

###
# Perform sum by ignoring NA
###
sumNoNA<-function(x){
  sum(x[which(!is.na(x))])
}

###
# Perform a two-class test of significance
###

two.class.test<-function(obs.table, all.exp, cond.control="C", cond.treatment="PA",to.plot=T){
  clus_names<-dimnames(obs.table)[[2]] #extract cluster names
  pp<-array(-1,length(clus_names)); 
  names(pp)<-clus_names;

  if(to.plot){
     par(mfrow=c(3,ceiling(length(clus_names)/3)));
  }
  for(this_clus in clus_names){
    #calculate difference in percentage in observed table
    obs.diff<-obs.table[cond.treatment,this_clus]/sum(obs.table[cond.treatment,]) - 
              obs.table[cond.control,this_clus]/sum(obs.table[cond.control,]); 
    #calculate difference in percentage for all permutated tables, list of 10000 differences in percentages
    all.diff<-all.exp[cond.treatment,this_clus,]/apply(all.exp[cond.treatment,,],2,sumNoNA) - 
              all.exp[cond.control,this_clus,]/apply(all.exp[cond.control,,],2,sumNoNA);
    if(to.plot){
       hist(all.diff,breaks=50,col="grey",border="grey",main=this_clus)
       abline(v=obs.diff,col="red")
    }
    
    # which is the index of the p that are different - calculates the lenght of the indexes = total number of percentages where observed is    bigger than permutated, do that for both bigger than permuted and smaller and permuted. Get something like 30% of permuted are smaller than observed 
    p.r<-length(which(obs.diff>all.diff))/length(which(!is.na(all.diff))); 
    p.l<-length(which(obs.diff<all.diff))/length(which(!is.na(all.diff)));
    pp[this_clus]<-min(p.r,p.l);
  }
  return(pp);
}