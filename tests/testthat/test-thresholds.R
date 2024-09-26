test_that("reparThr() output",{
  nthr <- 5

  # Michela fun
  mic_fun <-function(param)
  {
    ncat<-length(param)
    tmp<-param[-ncat]
    tmp[-1]<-exp(tmp[-1])
    param[-ncat]<-cumsum(tmp)
    param
  }

  set.seed(123)
  for (trial in 1:3) {
    param <- rnorm(nthr)
    thr <- reparThr(param, CON2UN = F)

    # test order
    expect_true(sum(sort(thr, decreasing = F) == thr)==length(thr))

    # test with Michela output
    expect_equal(thr, mic_fun(c(param, 1))[1:nthr])

    # check back and forth
    expect_equal(reparThr(thr, CON2UN = T), param)

  }
})
