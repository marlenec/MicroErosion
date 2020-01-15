##############################################################################################
#
# Function to calculate the loss of phylogenetic or taxonomic diversity of 
# a metacommunity (gamma-diversity) due to the disparition of local patches, based 
# on a user-supplied scenario
# e.g. as in Chiarello et al. ProcB, 2020, calculating expected loss of total microbiome diversity
# due to animal host species extinction.
#
# Function written by Marlène Chiarello (marlene.chiarello@gmail.com)
#
##############################################################################################

MicroErosion<-function(metacom=NULL, q=0, scenario=NULL, tree=NULL, r=1, nb_eq=F, example=F)
{
  require("entropart")
  require("ape")
  
  if(example==TRUE) 
	{
        if(!is.null(metacom)) warning("Erosion has been performed on example data. Your metacommunity has been ignored.")
	scenario<-rnorm(n = 10, mean = 1, sd = 0.2)
	metacom<-matrix(NA, nrow = 10, ncol = 10, dimnames=list(1:10, LETTERS[1:10]))
	names(scenario)<-rownames(metacom)
	for(i in 1:10)
	{
  		metacom[i, ]<-rnorm(n = 10, mean = 5, sd = 0.2)
	}
	tree<-rtree(n = 10, tip.label = LETTERS[1:10])
  } # end of example
  
  # check metacommunity matrix
  if(missing(metacom)) stop("A matrix with abundance data in metacommunity must be provided")
  metacom<-as.matrix(metacom)
  if(any(is.na(metacom))) stop("The abundances matrix contains NA(s). Please check")
  if(any(metacom<0)) stop("The abundances matrix contains negative values. Please check")
  metacom<-metacom[,which(colSums(metacom)>0), drop=F]

  # check replicates
    if(missing(r)) stop("A number of replicates must be provided")
    if(!is.numeric(r)) stop("The number of replicates is not numeric. Please check.")

  # Prepare vector of probabilities
  patches<-rownames(metacom)
  probas<-rep(1, length(patches))
  names(probas)<-patches
  if(!is.null(scenario))
  {
    # check scenario
    if(length(intersect(names(scenario), patches))<length(patches)) stop("Patches names in scenario and in metacommunity don't match. Please correct.")
    probas<-scenario
  } # end of if
  
  # Check tree
  if(!is.null(tree))
  {
    if(length(intersect(tree$tip.label, colnames(metacom)))<ncol(metacom)) stop("Species in tree and in metacommunity don't match. Please correct.")
  }

  # Preparing table to store results of Erosion scenario
  res<-matrix(data = NA, nrow=length(probas)-1, ncol=r, 
                dimnames=list(paste0((length(probas)-1):1, "_remaining"), 1:r))
  patches_removed<-res
  
  ###### Calculate initial diversity (same for all replicates)
  total_mat_rel<-colSums(metacom)/sum(colSums(metacom))
  print("Calculating initial diversity...")
	# Diversity
	if(is.null(tree))
      {
        div_ini<-Tsallis(total_mat_rel, q = q, CheckArguments=T) # Tsallis diversity of order q
        if (nb_eq==T) # transforming diversity into nb of eq species
        {
        	if(q==1) {div_ini<-exp(div_ini)}
        	if(q==2) {div_ini<-1/(1-div_ini)}
        } # end of if nbeq
  	} else
      {
          # prune tree
          todrop<-setdiff(tree$tip.label, names(total_mat_rel))
          tree<-drop.tip(phy = tree, tip = todrop, trim.internal = T)
          
          # Initial phylogenetic diversity
          div_ini<-ChaoPD(total_mat_rel, q=q, tree, CheckArguments=T, Normalize = nb_eq)
      } # end of else
  	names(div_ini)<-paste0("q_", q)
    print("...done. Patch erosion...")

  ###### Loop on replicates
  for(j in 1:r)
  {
    # Define proba, metacom, tree here before using them
    probas_j<-probas
    mat_j<-metacom
    if(!is.null(tree)){tree_j<-tree}
    
    # Iteration of patch loss
    for(i in 1:(length(probas)-1))
    {
      # Removal of 1 patch at a time in metacommunity table, according to its probability to get extinct
      patch_removed_i<-as.character(sample(x = names(probas_j), size = 1, prob = probas_j))
      probas_j<-probas_j[!(names(probas_j) == patch_removed_i)]
      patches_removed[i,j]<-as.character(patch_removed_i)
      mat_j<-mat_j[rownames(mat_j) %in% names(probas_j),, drop=F]
      mat_j<-mat_j[, which(colSums(mat_j) > 0), drop=F]
    
      # Get regional relative abundances
      total_mat<-colSums(mat_j)
      total_mat_rel<-total_mat/sum(total_mat)
    
      # Taxonomic diversity
      if(is.null(tree))
      {
        div_i<-Tsallis(total_mat_rel, q = q, CheckArguments=T) # Tsallis diversity of order q
        if (nb_eq==T) # transforming diversity into nb of eq species
        {
        	if(q==1) {div_i<-exp(div_i)}
        	if(q==2) {div_i<-1/(1-div_i)}
        } # end of if nbeq
  	  } # end of Taxonomic diversity
  	
  	# Phylogenetic diversity
  	else
  	{
  		# remove species from tree if needed
    	loss_i<-setdiff(names(total_mat_rel), colnames(metacom))
        tree_j<-drop.tip(phy = tree_j, tip = loss_i, trim.internal = T)
  	
  	 	div_i<-ChaoPD(total_mat_rel, q=q, tree_j, CheckArguments=T, Normalize = nb_eq)
  	} # end of Phylogenetic diversity
      
    # fill res table
    res[i,j]<-div_i
    } # end of i
    
    print(paste0(j, " out of ", r, " replicates done;"))
    
	} # end of j
      
  # storing values
  res_all<-list(div_ini, res, patches_removed)
  names(res_all)<-c("Initial_Gamma_diversity", "Gamma_diversity", "Patches_Removed")
  return(res_all)
	
} # end of MicroErosion

# To run example :
#q0<-MicroErosion(q = 0, r = 5, example=T)
#q1<-MicroErosion(q = 1, r = 5, nb_eq=F, example=T) # Allen's index of entropy. PDq as in Chao et al. 2010.
#q2<-MicroErosion(q = 1, r = 5, nb_eq=T, example=T) # Dq with q=2 as in Chao et al. 2010.

