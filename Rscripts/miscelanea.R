
# --- 1, 2, 3 sigmas

w <- qnorm(0.025)
1-2*pnorm(w)

plot_central_pnorm(q=w)

v <- qnorm(0.15)
1-2*pnorm(v)

plot_central_pnorm(q=v)

a <- qnorm(0.0015)
plot_central_pnorm(q=a)


# --- animate 

PATH <- paste0( getwd(), "/data/mohinora" )

listFILES <- list.files(path=PATH,
                        pattern = ".tif",
                        full.names = TRUE)

STACK <- stack(listFILES[c(2:4,1)])

reclassify(STACK, c(-Inf, -32768, NA, 32767, Inf, NA))

STACK <- STACK * 1e-4

animate(x=STACK,n=1)
