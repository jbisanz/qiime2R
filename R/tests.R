#tests.R

for(test in list.files("testfiles/", full.names = T)){
  message(test)
  tab<-read_qza(test)
}
