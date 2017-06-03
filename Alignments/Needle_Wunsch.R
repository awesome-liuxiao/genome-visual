needle.wunsch=function(sequence1,sequence2,gap,match,mismatch,print.maxscore=FALSE, print.matrix=FALSE)
{
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
  
  Alignment1=""
  Alignment2=""
  
  i=length(sequence1)
  j=length(sequence2)
  while(i>1 || j>1)
  {
    if (i > 1 && j > 1 &&  scoreMatrix[i,j] == (scoreMatrix[i-1,j-1] + mm(i,j) )    )
    {
      Alignment1 <- paste(sequence1[i],Alignment1,sep="")
      Alignment2 <- paste(sequence2[j],Alignment2,sep="")
      i = i - 1
      j = j - 1
    }
    else if ( (i > 1) && (scoreMatrix[i,j] == (scoreMatrix[i-1,j] + gap ) ) )
    {
      Alignment1 <- paste(sequence1[i],Alignment1,sep="")
      Alignment2 <- paste("-",Alignment2,sep="")
      i = i - 1
    }
    else if ( (j > 1)  && ( scoreMatrix[i,j] == (scoreMatrix[i,j-1] + gap) ) )
    {
      Alignment1 <- paste("-",Alignment1,sep="")
      Alignment2 <- paste(sequence2[j], Alignment2,sep="")
      j = j - 1
    }
  }
  
  if(print.matrix && print.maxscore)
  {
    output= list("Alignment1"=Alignment1,"Alignment2"= Alignment2, 
                 "maxscore"=max(scoreMatrix),"scoreMatrix"=scoreMatrix)
  }
  else if(!print.matrix && print.maxscore)
  {
    output=list("Alignment1"=Alignment1,"Alignment2"= Alignment2, 
                "maxscore"=max(scoreMatrix) ) 
  }
  else if(print.matrix && !print.maxscore)
  {
    output=list("Alignment1"=Alignment1,"Alignment2"= Alignment2, "scoreMatrix"=scoreMatrix) 
  }  
  else if(!print.matrix && !print.maxscore)
  {
    output=list("Alignment1"=Alignment1,"Alignment2"= Alignment2) 
  }
  return(output)
}

stringa= c("a test for the needleman-wunsch algorithm")
stringb= c("another test string for the needleman-wunsch algorithm to match")

needle.wunsch(stringa, stringb,match=1,mismatch=-1,gap=-1)
needle.wunsch(stringa,stringb,match=1,mismatch=0,gap=-1)

mousemyo= "MGLSDGEWQLVLNVWGKVEADLAGHGQEVLIGLFKTHPETLDKFDKFKNLKSEEDMKGSEDLKKHGCTVLTALGTILKKKGQHAAEIQPLAQSHATKHKIPVKYLEFISEIIIEVLKKRHSGDFGADAQGAMSKALELFRNDIAAKYKELGFQG"
humanmyo= "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"

needle.wunsch(mousemyo, humanmyo,match=1,mismatch=-1,gap=-1)
needle.wunsch(mousemyo,humanmyo,match=1,mismatch=0,gap=-1)