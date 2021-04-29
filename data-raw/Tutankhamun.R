## code to prepare `Tutankhamun` dataset goes here

library(pedtools)
Tutankhamun = readPed("data-raw/tutankhamun-jama.txt")

usethis::use_data(Tutankhamun, overwrite = TRUE)

#### Testing

showInTriangle(forrel::ibdEstimate(Tutankhamun)[11:6,], labels = T)

t1 = reconstruct(Tutankhamun, ids = labels(Tutankhamun)[1:4], noChildren = labels(Tutankhamun)[3:4], linear = FALSE, extra = 2)
plot(t1, 6)
# Try NorwFreq
q = 0.2; p = 1 - q
db = forrel::NorwegianFrequencies
db$D13S317 = c(p * db$D13S317, "16" = q)
db$D18S51 = c("8" = q, p * db$D18S51)
db$CSF1PO = c("6" = q, p * db$CSF1PO)
db$FGA = c(p * db$FGA, "31" = q)

x = setFreqDatabase(Tutankhamun, db)

showInTriangle(forrel::ibdEstimate(x)[11:6,], labels = T)

t1 = reconstruct(Tutankhamun, ids = labels(Tutankhamun)[3:6], linear = FALSE, extra = 2)


### AIC
ids = labels(Tutankhamun)[c(5,7)]
y = getAlleles(Tutankhamun, ids = ids)
ibdEstimate(Tutankhamun, ids)
