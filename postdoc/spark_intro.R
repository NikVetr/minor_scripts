Sys.getenv("JAVA_HOME")
Sys.setenv(JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home")
library(sparklyr)
# spark_available_versions()
# spark_install("2.0.0")
# spark_installed_versions()
spark_disconnect(sc)
sc <- spark_connect(master = "local", version = "2.0.0")
cars <- copy_to(sc, mtcars)
cars

# spark_web(sc)
library(DBI)
library(dplyr)
count(cars)
dbGetQuery(sc, "SELECT count(*) FROM mtcars")
select(cars, hp, mpg) %>%
  sample_n(100) %>%
  collect() %>%
  plot()
model <- ml_linear_regression(cars, mpg ~ hp)
model
model %>%
  ml_predict(copy_to(sc, data.frame(hp = 250 + 10 * 1:10))) %>%
  transmute(hp = hp, mpg = prediction) %>%
  full_join(select(cars, hp, mpg)) %>%
  collect() %>%
  plot()
