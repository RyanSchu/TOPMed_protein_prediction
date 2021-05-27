###data match sort
##takes in two data frames
##first checks if they are sorted the same, If yes, then done
##if not, by default will reorder input of --data to match input of --match
##assumes that the files contain the same samples, just in different orders

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparse)))

parser <- ArgumentParser()
parser$add_argument("--data", help="the data you would like to sort")
parser$add_argument("--match", help="the data set to match the order to")
parser$add_argument("--skip.col.data", help="number of columns of row names", type="integer")
parser$add_argument("--skip.col.match", help="number of columns of row names", type="integer")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

is.match<-function(names1,names2){
  if(length(names1) != length(names2)){
    cat("SampleSizeError: input data contain different number of samples\nCurrently not equipped to handle\n")
    quit("no",1)
  }
  else {
    compare<-sum(names1 == names2)
    check<-(compare == length(names1))
    return(check)
  }
}

check.fread<-function(file_name){
  zipped<-grepl(".gz$",file_name)
  if(zipped == T){
    file<-fread('zcat ' %&% file_name, header = T)
    return(file)
  } else{
    file<-fread(file_name, header = T)
    return(file)
  }
}

grab.ids<-function(file_name, start=1){
  df<-read.table(file_name,nrows = 1,header = F, stringsAsFactors = F)
  header<-unname(unlist(df[1,]))
  ids<-header[start:length(header)]
  return(ids)
}

data_ids<-grab.ids(args$data, (args$skip.col.data + 1))
match_ids<-grab.ids(args$match, (args$skip.col.match + 1))

if(is.match(data_ids,match_ids)){
  cat("file ids already in matching order\nExiting\n")
  quit("no",status=0)
} else {
  data<-check.fread(args$data)
  old_header<-colnames(data)
  final_header<-c(old_header[1:(args$skip.col.data)],match_ids)
  matched_data<-select(data,one_of(final_header))
  
  fwrite(matched_data,args$out,sep=' ',col.names = T,row.names = F,quote = F)
  cat(args$data, " reordered. Output as ", args$out, "\n")
}


#?quit
