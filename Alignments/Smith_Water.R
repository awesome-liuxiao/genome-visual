smith.water = function(sequence1,sequence2,gap,match,mismatch,print.maxscore=FALSE, print.matrix=FALSE) {
  
  if(length(strsplit(sequence1,"")[[1]])>1 && length(strsplit(sequence2,"")[[1]])>1) 
  {
    sequence1= strsplit(x=sequence1, split="")[[1]]
    sequence2= strsplit(x=sequence2, split="")[[1]]
  }

  sequence1=c("*",sequence1)
  sequence2=c("*",sequence2)
  
  mm=function(i,j)
  {
    if (sequence1[i]==sequence2[j])
      return(match)
    else
      return(mismatch)
  }
  
  costFunction = function(i,j)    
  {
    if(is.na(scoreMatrix[i,j])) 
    {
      if(i==j && i==1)            
      {
        return( 0 )
      }
      
      if(i==1)            
      {
        return( costFunction(i,j-1) + gap )
      }
      
      if(j==1)            
      {
        return( costFunction(i-1,j) + gap )
      }
      else            
      {
        costVector =c( ( costFunction(i-1,j) + gap ), 
                       ( costFunction(i,j-1) + gap ), 
                       ( costFunction(i-1, j-1) + mm(i,j) ) )
        return( costVector[which.max(costVector)] )
      }
    }
    else        
    {
      return (scoreMatrix[i,j])
    }
  }
  
  noRows=length(sequence1)
  noCols=length(sequence2)
  scoreMatrix=matrix(data=NA, nrow = noRows, ncol=noCols)
  rownames(scoreMatrix) = sequence1
  colnames(scoreMatrix) = sequence2
  
  for(i in 1:noRows)        
  {
    for (j in 1:noCols)        
    {
      scoreMatrix[i,j]=costFunction(i,j)
    }
  }
  
  max = which(scoreMatrix == max(scoreMatrix), arr.ind = TRUE)
  
  if(nrow(max) == 1) 
  {
    i=max[,1]
    j=max[,2]
  } 
  else 
  {
    for(k in 1:nrow(max)) 
    {
      i = max[,1]
      j = max[,2]
    }
  }
  
  Alignments1=c()
  Alignments2=c()
  
  for(n in length(i)) 
  {
    Alignment1=""
    Alignment2=""
    while(scoreMatrix[i[n] ,j[n] ] > 0)    
    {
      if (i[n] > 1 && j[n] > 1 &&  scoreMatrix[i[n],j[n]] == (scoreMatrix[i[n]-1,j[n]-1] + mm(i[n],j[n]) ))        
      {
        Alignment1 <- paste(sequence1[i[n]],Alignment1,sep="")
        Alignment2 <- paste(sequence2[j[n]],Alignment2,sep="")
        i[n] = i[n] - 1
        j[n] = j[n] - 1
      }
      else if ( (i[n] > 1) && (scoreMatrix[i[n],j[n]] == (scoreMatrix[i[n]-1,j[n]] + gap ) ) )        
      {
        Alignment1 <- paste(sequence1[i[n]],Alignment1,sep="")
        Alignment2 <- paste("-",Alignment2,sep="")
        i[n] = i[n] - 1
      }
      else if ( (j[n] > 1)  && ( scoreMatrix[i[n],j[n]] == (scoreMatrix[i[n],j[n]-1] + gap) ) )        
      {
        Alignment1 <- paste("-",Alignment1,sep="")
        Alignment2 <- paste(sequence2[j[n]], Alignment2,sep="")
        j[n] = j[n] - 1
      }
    }
    Alignments1 = c(Alignments1, Alignment1) 
    Alignments2 = c(Alignments2, Alignment2) 
  }
  if(print.matrix && print.maxscore)          
  {           
    output= list("Alignments1"=Alignments1,"Alignments2"= Alignments2, 
                 "maxscore"=max(scoreMatrix),"scoreMatrix"=scoreMatrix)
  }
  else if(!print.matrix && print.maxscore)   
  {
    output=list("Alignments1"=Alignments1,"Alignments2"= Alignments2, 
                "maxscore"=max(scoreMatrix) ) 
  }
  else if(print.matrix && !print.maxscore)    
  {
    output=list("Alignments1"=Alignments1,"Alignments2"= Alignments2, "scoreMatrix"=scoreMatrix) 
  }  
  else if(!print.matrix && !print.maxscore)    
  {
    output=list("Alignments1"=Alignments1,"Alignments2"= Alignments2) 
  }
  return(output)
}

stringa= c("a test for the needleman-wunsch algorithm")
stringb= c("another test string for the needleman-wunsch algorithm to match")

smith.water(stringa, stringb,match=1,mismatch=-1,gap=-1)
smith.water(stringa,stringb,match=1,mismatch=0,gap=-1)

mousemyo= "MGLSDGEWQLVLNVWGKVEADLAGHGQEVLIGLFKTHPETLDKFDKFKNLKSEEDMKGSEDLKKHGCTVLTALGTILKKKGQHAAEIQPLAQSHATKHKIPVKYLEFISEIIIEVLKKRHSGDFGADAQGAMSKALELFRNDIAAKYKELGFQG"
humanmyo= "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"

smith.water(mousemyo, humanmyo,match=1,mismatch=-1,gap=-1)
smith.water(mousemyo,humanmyo,match=1,mismatch=0,gap=-1)
