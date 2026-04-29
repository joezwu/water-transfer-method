args <- commandArgs(trailingOnly = F)
args

jobname <- sub("-","",args[length(args)])
jobname

datafile = sprintf("%s.outw",jobname)
data <- read.table(datafile)
nw <- data$V1
mean(nw)
sd(nw)/sqrt(length(nw))

pn <- as.data.frame(table( nw ))
ni <- as.numeric(as.matrix(pn)[,1])
N <- sum(pn$Freq)
pn
N
ni


q <- (pn$Freq/N)*choose(ni,1)
qm <- sum(q)
s <- (pn$Freq/N)*(choose(ni,1)-qm)^2
result <- sprintf("<n>            = %f +- %f", qm, sqrt(sum(s)/N))
write(result, "")

q <- (pn$Freq/N)*choose(ni,2)
qm <- sum(q)
s <- (pn$Freq/N)*(choose(ni,2)-qm)^2
result <- sprintf("<n_choose_2> = %f +- %f", qm, sqrt(sum(s)/N))
write(result, "")

q <- (pn$Freq/N)*choose(ni,3)
qm <- sum(q)
s <- (pn$Freq/N)*(choose(ni,3)-qm)^2
result <- sprintf("<n_choose_3> = %f +- %f", qm, sqrt(sum(s)/N))
write(result, "")

q <- (pn$Freq/N)*choose(ni,4)
qm <- sum(q)
s <- (pn$Freq/N)*(choose(ni,4)-qm)^2
result <- sprintf("<n_choose_4> = %f +- %f", qm, sqrt(sum(s)/N))
write(result, "")
