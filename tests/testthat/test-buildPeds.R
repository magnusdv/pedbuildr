bp = function(ids, sex = rep(1, length(ids)), verb = TRUE, ...) {
  buildPeds(ids, sex, verbose = verb, ...)
}

test_that("buildPeds cathces errors", {
  expect_error(bp(1:3, sex = 1:2, verb = F),
               "`ids` and `sex` must have the same length")
  expect_error(bp(1:3, sex = 1:3, verb = F),
               "Illegal elements found in `sex`")
  expect_error(bp(1:2, age = "a>b", verb = F),
               "Unknown ID label in `age`")
  expect_error(bp(1:2, age = "1", verb = F),
               "Character '>' missing in `age` entry")
})

test_that("buildPeds() parses age correctly", {
  expect_output(bp(1:3, age = c(1,1,1)), "Age info: -")
  expect_output(bp(1:3, age = c(1,1,1)), "Known non-PO: -", fixed = T)

  expect_output(bp(1:3, age = c(NA,NA,NA)), "Age info: -")
  expect_output(bp(1:3, age = c(NA,NA,NA)), "Known non-PO: -")

  expect_output(bp(1:3, age = c(2,1,1)), "Age info: 1>2, 1>3")
  expect_output(bp(1:3, age = c(2,NA,1)), "Age info: 1>3")
  expect_output(bp(1:3, age = c(2,2,1)), "Age info: 1>3, 2>3")

  expect_output(bp(1:3, age = "1>2,3"), "Age info: 1>2, 1>3")

})
