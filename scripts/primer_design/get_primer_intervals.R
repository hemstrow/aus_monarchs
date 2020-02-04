# this will create a file containing the primer target regions for primer3 SEQUEaCE_PRIMER_PAIR_OK_REGIOa_LIST argument.
e <- 7881485 # end of target region
s <- 7820867 # start of target region
o <- 500 # margin in which primers can be found
overlap <- 1500 # distance by which reads must overlap
sep <- 10000 # distance between each binding region, max allowable length of reads.
outfile <- "primer_regions.txt" # name of outfile





L <- e - s
L <- L - o*2
c <- sep
u <- overlap
max.spare <- 200
start <- 6500

# eq. 1: L + spare = ax - (u(a-1))
# eq. 2: x + 2o <= c
# eq. 3: spare <= max.spare
# x is the length per segment. a is the number of segments.
# eq 1 rearrange: a = ((L + spare) -u)/(x - u)
# eq2 rearrange: x <= c - 2o
# solve for a maximizing x, a must be an integer.
# spare can be any number, but should be minimzed (less than max.spare)
x <- c - 2*o
start <- x/2

opt.fun <- function(x, spare){
  return((((L + spare) - u))/(x - u))
}

# try every relevent option from start to x and from 1:max.spare
out <- outer(start:x, 1:max.spare, FUN = opt.fun)
colnames(out) <- 1:max.spare
rownames(out) <- start:x
out <- reshape2::melt(out)
out <- out[which(out[,3] %% 1 == 0),]
out <- out[which(out[,3] == min(out[,3])),]
best <- unlist(out[which.min(out[,2]),])
x <- best[1]
a <- best[3]

# get the sequence
st <- numeric(a)
st[1] <- 1
en <- numeric(a)
en[1] <- o + x

out <- character(a)
out[1] <- paste0(c(st[1], o, en[1], o), collapse = ",")

for(i in 2:(a - 1)){
  st[i] <- st[i - 1] + x - u
  en[i] <- en[i - 1]  - u + x
  out[i] <- paste0(c(st[i], o, en[i], o), collapse = ",")
}

st[a] <- st[a - 1] + x - u
en[a] <- (L + o*2) - o
out[a] <- paste0(c(st[a], o, en[a], o), collapse = ",")

writeLines(out, outfile)


