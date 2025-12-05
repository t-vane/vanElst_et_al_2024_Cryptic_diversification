### Script by Jordi Salmona

####EXAMPLE of "list_pairs" for a given SPECIES PAIR comparison "Mdan-Mdan" "Mrav-Mrav" "Mdan-Mrav"####


LOGnrme_regression<-function(df,list_pairs){
  
  #Split Pair Species Name into single Species name#
  pairs_sp = do.call(rbind,strsplit(list_pairs,'-'))
  
  #find WITHIN species comparisons#
  coord_with = which(pairs_sp[,1]==pairs_sp[,2])
  
  #find BETWEEN species comparisons#
  coord_betw = which(!pairs_sp[,1]==pairs_sp[,2])
  
  final_within = NULL
  
  ###Loop over each WITHIN species, separately###
  for(s in 1:length(coord_with)){
    
    sp_within_id = list_pairs[coord_with[s]]
    
    #subset WITHIN-Species dataset#
    sp_within_df = subset(df,combo==sp_within_id)
    
    ##perform linear regression between log-geo and genetic distance#
    lm.within = lm(value ~ log_geo_dist,data=sp_within_df)
    
    #get BETWEEN-Species ID#
    sp_between_id = list_pairs[coord_betw]
    #subset BEETWEEN-Species dataset#
    sp_between_df = subset(df,combo==sp_between_id)
    
    ###Predict the expected genetic distance given the BETWEEN-Species geographic distances.###
    ###The prediction is obtained from the "lm.within" WITHIN-Species linear regression parameters##
    yy<-as.vector(predict.lm(lm.within, sp_between_df))
    
    ###It computes NRMSE using deviation of the OBSERVED VALUES ("sp_between_df$value") from###
    ###the predicted values ("yy"). Normalized by the range of observed values###
    nrmse<-sqrt(mean((yy - sp_between_df$value)^2,na.rm=TRUE)) / diff(range(sp_between_df$value,na.rm = TRUE))
    
    tmp_within<-data.frame(within=sp_within_id,between=sp_between_id,nrmse=nrmse)
    final_within<-rbind(final_within,tmp_within)
  }
  
  return(final_within)
}
