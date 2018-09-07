censdat <- data.frame(read.csv("sc-est2016-alldata5.csv", header = T))

censdat <- censdat[which(censdat$SEX==0 & censdat$ORIGIN==0),]

st_list <- list(c('Connecticut','Maine','Massachusetts','New Hampshire','Rhode Island','Vermont'),
                c('New Jersey','New York','Puerto Rico','The Virgin Islands'), 
                c('Delaware','District of Columbia','Maryland','Pennsylvania','Virginia',
                  'West Virginia'),
                c('Alabama','Florida','Georgia','Kentucky','Mississippi','North Carolina',
                  'South Carolina','Tennessee'),
                c('Illinois','Indiana','Michigan','Minnesota','Ohio','Wisconsin'),
                c('Arkansas','Louisiana','New Mexico','Oklahoma','Texas'),
                c('Iowa','Kansas','Missouri','Nebraska'),
                c('Colorado','Montana','North Dakota','South Dakota','Utah','Wyoming'),
                c('Arizona','California','Hawaii','Nevada','American Samoa',
                  'Commonwealth of the Northern Mariana Islands','Federated States of Micronesia',
                  'Guam','Marshall Islands','Republic of Palau'),
                c('Alaska','Idaho','Oregon','Washington'))


age_list <- list(c(0:4),c(5:19),c(20:64),c(65:85))
HHSagepop <- NULL
for (HHSreg in 1:10){
  agepop <- NULL
  regselind <- which(censdat$NAME %in% st_list[[HHSreg]])
  tempdat <- censdat[regselind,]
  for (ag in 1:4){
    agselind <- which(tempdat$AGE %in% age_list[[ag]])
    agepop <- c(agepop,sum(tempdat$POPESTIMATE2016[agselind]))
  }
  HHSagepop <- rbind(HHSagepop,c(HHSreg, agepop),deparse.level = 0)
}

outfilename <- paste0('pop_HHS_4.dat')
write.table(HHSagepop,outfilename,quote = F,col.names = F,row.names = F)

