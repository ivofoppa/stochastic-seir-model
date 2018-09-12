censdat <- data.frame(read.csv("sc-est2016-alldata5.csv", header = T))

censdat <- censdat[which(censdat$SEX==0 & censdat$ORIGIN==0),]

st_list <- c('Connecticut','Maine','Massachusetts','New Hampshire','Rhode Island','Vermont',
                'New Jersey','New York','Puerto Rico','The Virgin Islands', 
                'Delaware','District of Columbia','Maryland','Pennsylvania','Virginia',
                  'West Virginia',
                'Alabama','Florida','Georgia','Kentucky','Mississippi','North Carolina',
                  'South Carolina','Tennessee',
                'Illinois','Indiana','Michigan','Minnesota','Ohio','Wisconsin',
                'Arkansas','Louisiana','New Mexico','Oklahoma','Texas',
                'Iowa','Kansas','Missouri','Nebraska',
                'Colorado','Montana','North Dakota','South Dakota','Utah','Wyoming',
                'Arizona','California','Hawaii','Nevada','American Samoa',
                  'Commonwealth of the Northern Mariana Islands','Federated States of Micronesia',
                  'Guam','Marshall Islands','Republic of Palau',
                'Alaska','Idaho','Oregon','Washington')


age_list <- list(c(0:4),c(5:19),c(20:64),c(65:85))
stateagepop <- NULL
for (state in st_list){
  stateselind <- which(censdat$NAME %in% state)
  if (length(stateselind) > 0 ){
    tempdat <- censdat[stateselind,]
    
    agepop <- NULL
    for (ag in 1:4){
      agselind <- which(tempdat$AGE %in% age_list[[ag]])
      agepop <- c(agepop,sum(tempdat$POPESTIMATE2016[agselind]))
    }
    stateagepop <- rbind(stateagepop,c(state, agepop),deparse.level = 0)
  }
}

stateagepop <- data.frame(stateagepop)
colnames(stateagepop) <- c('state','ag1','ag2','ag3','ag4')

for (ag in 1:4){
  eval(parse(text = paste0('stateagepop$ag',ag,' <- as.numeric(as.vector(stateagepop$ag',ag,'))')))
}

outfilename <- paste0('pop_state_4.csv')
write.csv(stateagepop,outfilename)



