library(relater)
library(ggplot2)


### load relate coal files and plot effective popsizes later


coal <- read.coal("relate_popsize_CAN_FR.pairwise_0kb.coal")

coal$popsize <- 0.5/coal$haploid.coalescence.rate
coal$years<-1.5*(coal$epoch.start)
coal$within <- c("across groups","within group")[(coal$group1 == coal$group2)+1]

pop_select <- c("ES03", "ES04", "ES17", "ES25", "FR01")
coal_select <- coal[coal$group1%in%pop_select & coal$group2%in%pop_select,]
ggplot(coal_select)+
  geom_step(aes(x = epoch.start, y = popsize, color = group2, linetype = within))+
  scale_x_continuous(trans = "log10", limit = c(1e3, 1e6)) +
  scale_y_continuous(trans = "log10", limit = c(1e3,1e6)) +
  annotation_logticks(sides = "lb") +
  scale_linetype_manual(values = c(2,1), name = "") +
  facet_grid(rows = .~group1) +
  xlab("years ago") +
  ylab("diploid effective population size")

ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_step(aes(x = years, y = popsize, color=group1))+
  scale_x_continuous(trans = "log10", limit = c(1e3, 1e6)) +
  scale_y_continuous(trans = "log10", limit = c(1e3,1e6)) +
  annotation_logticks(sides = "lb") +
  scale_linetype_manual(values = c(2,1), name = "") +
  xlab("years ago") +
  ylab("diploid effective population size")


library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

samples_sizes
unique(coal2$years)
ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_step(aes(x = years, y = round(popsize, digits = 2), color = group1))+
  scale_x_continuous(trans = reverselog_trans(10), limits = c(2682690, 1e3), breaks = c(1e6, 1e5,1e4,1e3), labels = c(expression(paste(10^6)),expression(paste(10^5)),expression(paste(10^4)),expression(paste(10^3))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  scale_linetype_manual(values = c(2,1), name = "") +
  xlab("Years ago") +
  ylab(expression(paste("Effective diploid population size [x", 10^3, "]",sep="")))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,2,0,1), "cm"))+
  labs(color="Population")


ggplot(coal_select[coal_select$group1%in%c("FR01", "FR08"),])+
  geom_step(aes(x = years, y = round(popsize, digits = 2), color = group2, linetype=within))+
  scale_x_continuous(trans = reverselog_trans(10), limits = c(2682690, 1e3), breaks = c(1e6, 1e5,1e4,1e3), labels = c(expression(paste(10^6)),expression(paste(10^5)),expression(paste(10^4)),expression(paste(10^3))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb") +
  scale_linetype_manual(values = c(2,1), name = "") +
  xlab("Years ago") +
  facet_grid(rows = .~group1) +
  ylab(expression(paste("Effective population size (x", 10^3, ")",sep="")))+
  theme(axis.text = element_text(size=10), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=17), legend.title = element_text(size=17), legend.key.size = unit(1, "cm"), legend.text = element_text(size = 15))+
  labs(color="Population")



#load bootstrap results
head(coal_boot)
coal_boot <- read.coal(paste0("bootstraps_0kb/relate_popsize_CAN_FR_boot",as.character(1),".pairwise.coal"))
# coal_boot$popsize <- 0.5/coal_boot$haploid.coalescence.rate
coal_boot$years<-1.5*(coal_boot$epoch.start)
coal_boot$within <- c("across groups","within group")[(coal_boot$group1 == coal_boot$group2)+1]
colnames(coal_boot)[4] <- "haploid.coalescence.rate.boot1"
coal_boot<- coal_boot[,c(1:3,5:6,4)]
coal_boot$rCCR<-NA
for(i in 1:nrow(coal_boot)){
  coal_boot$rCCR[i] <- 2*coal_boot$haploid.coalescence.rate.boot1[i]/(coal_boot$haploid.coalescence.rate.boot1[coal_boot$group1== coal_boot$group1[i]& coal_boot$group2==coal_boot$group1[i] & coal_boot$epoch.start==coal_boot$epoch.start[i]]+coal_boot$haploid.coalescence.rate.boot1[coal_boot$group1== coal_boot$group2[i]& coal_boot$group2==coal_boot$group2[i] & coal_boot$epoch.start==coal_boot$epoch.start[i]] )
}
colnames(coal_boot)[5] <- "within"
colnames(coal_boot)[7]<- "rCCR.boot1"
head(coal_boot)

for (i in 2:100){
  print(i)
  coal_tmp <- read.coal(paste0("bootstraps_0kb/relate_popsize_CAN_FR_boot",as.character(i),".pairwise.coal"))
  coal_tmp$rCCR<-NA
  for(j in 1:nrow(coal_boot)){
    coal_tmp$rCCR[j] <- 2*coal_tmp$haploid.coalescence.rate[j]/(coal_tmp$haploid.coalescence.rate[coal_tmp$group1== coal_tmp$group1[j]& coal_tmp$group2==coal_tmp$group1[j] & coal_tmp$epoch.start==coal_tmp$epoch.start[j]]+coal_tmp$haploid.coalescence.rate[coal_tmp$group1== coal_tmp$group2[j]& coal_tmp$group2==coal_tmp$group2[j] & coal_tmp$epoch.start==coal_tmp$epoch.start[j]] )
  }
  coal_boot <- cbind(coal_boot, coal_tmp$haploid.coalescence.rate,coal_tmp$rCCR)
  colnames(coal_boot)[c((ncol(coal_boot)-1), ncol(coal_boot))] <-c(paste0( "haploid.coalescence.rate.boot", as.character(i)),paste0( "rCCR.boot", as.character(i)))
}


head(coal_boot)

test2<-coal_boot[coal_boot$group1=="FR08" & coal_boot$group2=="FR08" ,]
test2<- test2[test2$years==100000.000,]
sort(test2[1,grep(pattern = "haploid.coalescence.rate.boot", x = colnames(test2))][!is.na(test2[1,grep(pattern = "haploid.coalescence.rate.boot", x = colnames(test2))])])
head(test)



coal_boot$haploid.coalescence.rate.min <- NA
coal_boot$haploid.coalescence.rate.max <- NA
coal_boot$rCCR.min <- NA
coal_boot$rCCR.max <- NA
colnames(coal_boot)=="within group"
i<-1
is.na(coal_boot$haploid.coalescence.rate.max)
coal_boot[coal_boot$group1=="ES01"& coal_boot$group2=="FR18",]
for(i in 1:nrow(coal_boot)){
  tmp<- sort(coal_boot[i,grep(pattern = "haploid.coalescence.rate.boot", x = colnames(coal_boot))][!is.na(coal_boot[i,grep(pattern = "haploid.coalescence.rate.boot", x = colnames(coal_boot))])])
  tmp2<- sort(coal_boot[i,grep(pattern = "rCCR.boot", x = colnames(coal_boot))][!is.na(coal_boot[i,grep(pattern = "rCCR.boot", x = colnames(coal_boot))])])
  if(length(tmp)>20){
    coal_boot$haploid.coalescence.rate.min[i] <- tmp[round(0.05*length(tmp))]
    coal_boot$haploid.coalescence.rate.max[i] <- tmp[round(0.95*length(tmp))]
  }
  if(length(tmp2)>20){
    coal_boot$rCCR.min[i] <- tmp2[round(0.05*length(tmp2))]
    coal_boot$rCCR.max[i] <- tmp2[round(0.95*length(tmp2))]
  }
}

#write.table(x = coal_boot, file = "coal_boot_0kb.txt", sep="\t", quote = F, col.names = T, row.names = F)
head(coal_boot)
coal_boot<-read.table("coal_boot_0kb.txt", header = T, sep="\t")
head(coal_boot)

coal$haploid.coalescence.rate.min<- coal_boot$haploid.coalescence.rate.min
coal$haploid.coalescence.rate.max<- coal_boot$haploid.coalescence.rate.max
coal$popsize_max <- 0.5/coal$haploid.coalescence.rate.min
coal$popsize_min <- 0.5/coal$haploid.coalescence.rate.max
coal$rCCR.min<- coal_boot$rCCR.min
coal$rCCR.max<- coal_boot$rCCR.max

coal$rCCR<-NA
for(j in 1:nrow(coal)){
  coal$rCCR[j] <- 2*coal$haploid.coalescence.rate[j]/(coal$haploid.coalescence.rate[coal$group1== coal$group1[j]& coal$group2==coal$group1[j] & coal$epoch.start==coal$epoch.start[j]]+coal$haploid.coalescence.rate[coal$group1== coal$group2[j]& coal$group2==coal$group2[j] & coal$epoch.start==coal$epoch.start[j]] )
}



####
coal_select <- coal[coal$group1%in%pop_select & coal$group2%in%pop_select,]
unique(coal_select$group1)
unique(coal_select$group2)

test<-coal_select[coal_select$group1=="FR08" &coal_select$group2=="FR08"  ,]
colnames(test)
test[,c(6,10,5,11)]

ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_ribbon(aes(ymin = popsize_min, ymax = popsize_max, x=years,fill = group1), alpha = 0.2) + # Shaded CI area
  geom_path(aes(x = years, y = round(popsize, digits = 2), color = group1),size=1.2)+
  scale_x_continuous(trans = reverselog_trans(10), limits = c(2682690, 1e3), breaks = c(1e6, 1e5,1e4,1e3), labels = c(expression(paste(10^6)),expression(paste(10^5)),expression(paste(10^4)),expression(paste(10^3))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3.5e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  scale_linetype_manual(values = c(2,1), name = "") +
  scale_color_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"), labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  scale_fill_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"),labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  xlab("Years ago") +
  ylab(expression(paste("Effective diploid population size [x", 10^3, "]",sep="")))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,2,0,1), "cm"))+
  labs(color="Population")+
  labs(color="Population", fill="Population")

stairstepn <- function( data, direction="hv", yvars="y" ) {
  direction <- match.arg( direction, c( "hv", "vh" ) )
  data <- as.data.frame( data )[ order( data$x ), ]
  n <- nrow( data )
  
  if ( direction == "vh" ) {
    xs <- rep( 1:n, each = 2 )[ -2 * n ]
    ys <- c( 1, rep( 2:n, each = 2 ) )
  } else {
    ys <- rep( 1:n, each = 2 )[ -2 * n ]
    xs <- c( 1, rep( 2:n, each = 2))
  }
  
  data.frame(
    x = data$x[ xs ]
    , data[ ys, yvars, drop=FALSE ]
    , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
  ) 
}

stat_stepribbon <- 
  function(mapping = NULL, data = NULL, geom = "ribbon", position = "identity", inherit.aes = TRUE) {
    ggplot2::layer(
      stat = Stepribbon, mapping = mapping, data = data, geom = geom, 
      position = position, inherit.aes = inherit.aes
    )
  }

StatStepribbon <- 
  ggproto("stepribbon", Stat,
          compute_group = function(., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
            stairstepn( data = data, direction = direction, yvars = yvars )
          },                        
          required_aes = c( "x", "ymin", "ymax" )
  )



ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_ribbon(aes(ymin = popsize_min, ymax = popsize_max, x=years,fill = group1), alpha = 0.2, stat = "stepribbon") + # Shaded CI area
  geom_step(aes(x = years , y = round(popsize, digits = 2), color = group1),size=1.2)+
  scale_x_continuous(trans = reverselog_trans(10), limits = c(2682690, 1e3), breaks = c(1e6, 1e5,1e4,1e3), labels = c(expression(paste(10^6)),expression(paste(10^5)),expression(paste(10^4)),expression(paste(10^3))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3.5e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  scale_linetype_manual(values = c(2,1), name = "") +
  scale_color_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"), labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  scale_fill_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"),labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  xlab("Years ago") +
  ylab(expression(paste("Effective diploid population size [x", 10^3, "]",sep="")))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,2,0,1), "cm"))+
  labs(color="Population", fill="Population")


##look backwards in time
max(coal_select$years)
popsize<-ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_ribbon(aes(ymin = popsize_min, ymax = popsize_max, x=years,fill = group1), alpha = 0.2, stat = "stepribbon") + # Shaded CI area
  geom_step(aes(x = years , y = round(popsize, digits = 2), color = group1),size=1.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3.5e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  scale_linetype_manual(values = c(2,1), name = "") +
  scale_color_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"), labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  scale_fill_manual(breaks=c("ES03","ES04", "ES17", "FR01", "FR08"),labels=c("ES03","ES04", "ES17", "FR01", "FR08"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  xlab("Years ago") +
  ylab(expression(paste("Effective population size [x", 10^3, "]",sep="")))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.15, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,2,0,1), "cm"))+
  labs(color="Population", fill="Population")+
  annotate(geom = "rect",xmin=c(1e4,19e3), xmax=c(1.15e5, 26.5e3), ymin=c(-Inf,-Inf), ymax=c(Inf,Inf), fill="blue", alpha=0.2)

popsize<-ggplot(coal_select[coal_select$group1==coal_select$group2,])+
  geom_ribbon(aes(ymin = popsize_min, ymax = popsize_max, x=years,fill = group1), alpha = 0.2, stat = "stepribbon") + # Shaded CI area
  geom_step(aes(x = years , y = round(popsize, digits = 2), color = group1),size=1.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  scale_y_continuous(limit = c(0, 3.5e5), breaks =seq(0,3e5, by=1e5), labels= seq(0,300, by=100)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  scale_linetype_manual(values = c(2,1), name = "") +
  scale_color_manual(breaks=c("ES03","ES04", "ES17","ES25", "FR01"), labels=c("ES03","ES04", "ES17", "ES25", "FR01"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  scale_fill_manual(breaks=c("ES03","ES04", "ES17","ES25", "FR01"),labels=c("ES03","ES04", "ES17", "ES25", "FR01"),values = c("green4", "#6A3D9A", "maroon", "#56B4E9", "#BCBD22"))+
  xlab("Years ago") +
  ylab(expression(paste("Effective population size [x", 10^3, "]",sep="")))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.15, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,2,0,1), "cm"))+
  labs(color="Population", fill="Population")+
  annotate(geom = "rect",xmin=c(1e4,19e3), xmax=c(1.15e5, 26.5e3), ymin=c(-Inf,-Inf), ymax=c(Inf,Inf), fill="blue", alpha=0.2)

#### Calculate rCCR values
rCCR_only<-c()
i<-1
for(i in 1:length(unique(coal$group1))){
  tmp<-coal[coal$group1==unique(coal$group1)[i],]
  tmp2<-tmp[tmp$group2==unique(coal$group1)[i],]
  for(j in 1:length(unique(tmp$group2))){
    if(unique(tmp$group2)[j]!=unique(coal$group1)[i]){
      tmp3<-tmp[tmp$group2==(unique(tmp$group2)[j]),]
      tmp4<-coal[coal$group1==unique(tmp$group2)[j] & coal$group2==unique(tmp$group2)[j],]
      tmp3$rCCR<-2*tmp3$haploid.coalescence.rate/(tmp2$haploid.coalescence.rate+ tmp4$haploid.coalescence.rate)
      rCCR_only<-rbind(rCCR_only,tmp3)
    }
  }
}

head(rCCR_only)
coal <- rCCR_only

sum(is.na(coal$region1))
sum(is.na(coal$region2))
coal[is.na(coal$region1),]
unique(coal$group1[is.na(coal$region1)])
unique(coal$group1[is.na(coal$region2)])



rCCR_tab<-c()
i<-1
j<-1
k<-1
for(i in 1:length(unique(coal$group1))){
  tmp<-coal[coal$group1==unique(coal$group1)[i],]
   tmp2<-tmp[tmp$group2==unique(coal$group1)[i],]
  for(j in 1:length(unique(tmp$group2))){
    if(unique(tmp$group2)[j]!=unique(coal$group1)[i]){
      tmp3<-tmp[tmp$group2==(unique(tmp$group2)[j]),]
      rCCR0.25<-"Na"
      tCCR0.25<-"Na"
      rCCR0.5<-"Na"
      tCCR0.5<-"Na"
      rCCR0.75<-"Na"
      tCCR0.75<-"Na"
      t1<-TRUE
      t2<-TRUE
      for(k in 1:nrow(tmp3)){
        if(is.na(tmp3$rCCR[k])){
          next;
        }
        if(tmp3$rCCR[k]>0.25 & t1 & k!=1){
          rCCR0.25<-tmp3$rCCR[(k-1)]
          tCCR0.25<-tmp3$years[(k-1)]
          t1<-FALSE
          if(tmp3$rCCR[k]>0.5 & t2){
            rCCR0.5<-tmp3$rCCR[(k-1)]
            tCCR0.5<-tmp3$years[(k-1)]
            t2<-FALSE
          }
        }else if(tmp3$rCCR[k]>0.5 & t2  & k!=1){
          rCCR0.5<-tmp3$rCCR[(k-1)]
          tCCR0.5<-tmp3$years[(k-1)]
          t2<-FALSE
          if(tmp3$rCCR[k]>0.75 & !is.na(tmp3$rCCR[k])){
            rCCR0.75<-tmp3$rCCR[(k-1)]
            tCCR0.75<-tmp3$years[(k-1)]
            tmp4<-c(unique(coal$group1)[i], coal$region1[coal$group1==unique(coal$group1)[i]][1],unique(tmp$group2)[j],tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.25, tCCR0.25, "0.25")
            rCCR_tab<-rbind(rCCR_tab, tmp4)
            tmp4<-c(unique(coal$group1)[i],coal$region1[coal$group1==unique(coal$group1)[i]][1],unique(tmp$group2)[j], tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.5, tCCR0.5,"0.5")
            rCCR_tab<-rbind(rCCR_tab, tmp4)
            tmp4<-c(unique(coal$group1)[i],coal$region1[coal$group1==unique(coal$group1)[i]][1], unique(tmp$group2)[j], tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.75,  tCCR0.75, "0.75")
            rCCR_tab<-rbind(rCCR_tab, tmp4)
            print(i)
            break  
          }
        }else if(tmp3$rCCR[k]>0.75 & !is.na(tmp3$rCCR[k])  & k!=1 | k==nrow(tmp3)){
          rCCR0.75<-tmp3$rCCR[(k-1)]
          tCCR0.75<-tmp3$years[(k-1)]
          tmp4<-c(unique(coal$group1)[i],coal$region1[coal$group1==unique(coal$group1)[i]][1], unique(tmp$group2)[j], tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.25, tCCR0.25,"0.25")
          rCCR_tab<-rbind(rCCR_tab, tmp4)
          tmp4<-c(unique(coal$group1)[i],coal$region1[coal$group1==unique(coal$group1)[i]][1],unique(tmp$group2)[j], tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.5, tCCR0.5, "0.5")
          rCCR_tab<-rbind(rCCR_tab, tmp4) 
          tmp4<-c(unique(coal$group1)[i],coal$region1[coal$group1==unique(coal$group1)[i]][1], unique(tmp$group2)[j], tmp$region2[tmp$group2==unique(tmp$group2)[j]][1], rCCR0.75,  tCCR0.75, "0.75")
          rCCR_tab<-rbind(rCCR_tab, tmp4)
          print(i)
          break
        }
      }
    }
  }
}


rCCR_tab<-as.data.frame(rCCR_tab)
colnames(rCCR_tab)<-c("Pop1","Pop2", "rCCR", "years", "cutoff")
rCCR_tab$years<-as.numeric(rCCR_tab$years)
rCCR_tab<-rCCR_tab[rCCR_tab$Pop1!=rCCR_tab$Pop2,]


write.table(rCCR_tab, file = "rCCR_tab.txt", sep="\t", col.names = T, row.names = F, quote = F)


###rCCR with bootstrap and local polynomial regression to infer split times at rCCR=0.5



sum(is.na(coal$region1))
sum(is.na(coal$region2))
coal[is.na(coal$region1),]
unique(coal$group1[is.na(coal$region1)])
unique(coal$group1[is.na(coal$region2)])

library(mgcv)

head(coal)
split_times<-c()
i<-1
j<-2
k<-2
coal <- coal[coal$group1!="ES15" & coal$group2!="ES15" & coal$group1!="FR18" & coal$group2!="FR18",]
for(i in 1:length(unique(coal$group1))){
  tmp<-coal[coal$group1==unique(coal$group1)[i],]
  #  tmp2<-tmp[tmp$group2==unique(coal$group1)[i],]
  for(j in 1:length(unique(tmp$group2))){
    if(unique(tmp$group2)[j]!=unique(coal$group1)[i]){
      tmp2<-tmp[tmp$group2==(unique(tmp$group2)[j]),]
      
      tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
      log_years<- log(tmp2$years)
      spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
      spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df = 20)  # Adjust df for more or less flexibility
      spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df = 20)  # Adjust df for more or less flexibility
      max_years<- max(log_years)
      min_years<- min(log_years)
      
      new_years<- seq(min_years,max_years, length.out=1000)
      predicted_rCCR <- predict(spline_fit,x=new_years)
      predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
      predicted_rCCR_max <- predict(spline_fit_max, x=new_years)
      
      # plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
      # lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
      # lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
      # lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

      closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
      closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
      closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))
      
      
      split_times<- rbind(split_times, c(unique(coal$group1)[i],unique(tmp$group2)[j],exp(predicted_rCCR$x[closest_index]),exp(predicted_rCCR$x[closest_index_min]),exp(predicted_rCCR$x[closest_index_max])))
    }
  }
}

##split times within groups 

##

split_times<- as.data.frame(split_times)
split_times$point_estimate <- as.numeric(split_times$point_estimate)
colnames(split_times)<- c("pop1", "pop2", "point_estimate", "max", "min")
# write.table(split_times, file = "Split_times_0kb.txt", sep="\t", col.names = T,row.names = F, quote=F)
tmp1<-split_times[grep(split_times$pop1, pattern = "FR"),]
split_times_FR <- tmp1[grep(tmp1$pop2, pattern = "FR"),]
split_times_FR$point_estimate <- as.numeric(split_times_FR$point_estimate)
max(split_times_FR$point_estimate, na.rm=T)

split_times[split_times$pop1=="FR01" & split_times$pop2=="FR02",]
split_times[split_times$pop1=="FR01" & split_times$pop2=="ES02",]
split_times[split_times$pop1=="FR02" & split_times$pop2=="ES04",]
split_times[split_times$pop1=="FR10" & split_times$pop2=="ES04",]
split_times[split_times$pop1=="ES17" & split_times$pop2=="ES25",]
split_times[split_times$pop1=="ES03" & split_times$pop2=="ES25",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="ES25",]

split_times[split_times$pop1=="ES17" & split_times$pop2=="ES25",]
split_times[split_times$pop1=="ES03" & split_times$pop2=="ES25",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="ES25",]

split_times[split_times$pop1=="ES17" & split_times$pop2=="ES24",]
split_times[split_times$pop1=="ES03" & split_times$pop2=="ES24",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="ES24",]

split_times[split_times$pop1=="ES04" & split_times$pop2=="ES24",]


split_times[split_times$pop1=="FR02" & split_times$pop2=="FR10",]


split_times[split_times$pop1=="FR02" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="FR07" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="FR10" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="FR19" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="FR19" & split_times$pop2=="FR10",]

split_times[split_times$pop1=="FR20" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="FR21" & split_times$pop2=="FR01",]

split_times[split_times$pop1=="FR20" & split_times$pop2=="FR19",]

split_times[split_times$pop1=="FR19" & split_times$pop2=="FR01",]


split_times_group1 <- split_times[split_times$pop1%in%c("ES23", "ES05", "ES04", "ES10") &split_times$pop2%in%c("ES23", "ES05", "ES04", "ES10"), ]
split_times_group1$point_estimate
min(split_times_group1$point_estimate)
max(split_times_group1$point_estimate)


split_times_group2 <- split_times[split_times$pop1%in%c("ES01", "ES02", "ES03", "ES06") &split_times$pop2%in%c("ES01", "ES02", "ES03", "ES06"), ]
max(split_times_group2$point_estimate)
min(split_times_group2$point_estimate)


split_times_group3 <- split_times[split_times$pop1%in%c("ES08", "ES09", "ES17", "ES24", "ES25") &split_times$pop2%in%c("ES08", "ES09", "ES17", "ES24", "ES25"), ]
split_times_group3$point_estimate
max(split_times_group3$point_estimate)
min(split_times_group3$point_estimate)
split_times[split_times$pop1=="ES03" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="ES17" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="ES24" & split_times$pop2=="FR01",]
split_times[split_times$pop1=="ES25" & split_times$pop2=="FR01",]

split_times[split_times$pop1=="ES03" & split_times$pop2=="FR08",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="FR08",]
split_times[split_times$pop1=="ES17" & split_times$pop2=="FR08",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="ES03",]
split_times[split_times$pop1=="ES04" & split_times$pop2=="ES17",]
split_times[split_times$pop1=="ES03" & split_times$pop2=="ES17",]
split_times$point_estimate<- as.numeric(split_times$point_estimate)

split_times_ES <- split_times[grep(split_times$pop1,pattern= "ES"),]
split_times_ES <- split_times_ES[grep(split_times_ES$pop2,pattern= "ES"),]
split_times_ES$point_estimate <- as.numeric(split_times_ES$point_estimate)
min(split_times_ES$point_estimate)
max(split_times_ES$point_estimate)
split_times_ES_W_to_CE<- split_times_ES[split_times_ES$pop1 %in% c("ES17","ES08","ES09", "ES24","ES25"),]
split_times_ES_W_to_CE<-split_times_ES_W_to_CE[!split_times_ES_W_to_CE$pop2 %in% c("ES17","ES08","ES09", "ES24","ES25"), ]
max(split_times_ES_W_to_CE$point_estimate)
min(split_times_ES_W_to_CE$point_estimate)

split_times_ES_W_to_W<- split_times_ES[!split_times_ES$pop1 %in% c("ES17","ES08","ES09", "ES24","ES25"),]
split_times_ES_W_to_W<-split_times_ES_W_to_W[!split_times_ES_W_to_W$pop2 %in% c("ES17","ES08","ES09", "ES24","ES25"), ]
split_times_ES_W_to_W<-split_times_ES_W_to_W[split_times_ES_W_to_W$pop1 %in% c("ES04","ES05","ES23"), ]
split_times_ES_W_to_W<-split_times_ES_W_to_W[!split_times_ES_W_to_W$pop2 %in% c("ES04","ES05","ES23"), ]

max(split_times_ES_W_to_W$point_estimate)
min(split_times_ES_W_to_W$point_estimate)

max(split_times_ES_W_to_CE$point_estimate)/min(split_times_ES_W_to_W$point_estimate)
min(split_times_ES_W_to_CE$point_estimate)/max(split_times_ES_W_to_W$point_estimate)
mean(split_times_ES_W_to_CE$point_estimate)/mean(split_times_ES_W_to_W$point_estimate)

split_times_ES_CE<- split_times_ES[split_times_ES$pop1 %in% c("ES17","ES08","ES09", "ES24","ES25"),]
split_times_ES_CE<- split_times_ES_CE[split_times_ES_CE$pop2 %in% c("ES17","ES08","ES09", "ES24","ES25"),]
min(split_times_ES_CE$point_estimate)
max(split_times_ES_CE$point_estimate)

split_times_ES_W2<- split_times_ES[split_times_ES$pop1 %in% c("ES04","ES05","ES23"),]
split_times_ES_W2<- split_times_ES_W2[split_times_ES_W2$pop2 %in% c("ES04","ES05","ES23"),]
min(split_times_ES_W2$point_estimate)
max(split_times_ES_W2$point_estimate)


split_times_ES_W1<- split_times_ES[split_times_ES$pop1 %in% c("ES01","ES02","ES03", "ES06", "ES10"),]
split_times_ES_W1<- split_times_ES_W1[split_times_ES_W1$pop2 %in% c("ES01","ES02","ES03", "ES06", "ES10"),]
min(split_times_ES_W1$point_estimate)
max(split_times_ES_W1$point_estimate)



split_times_ES[which(split_times_ES$point_estimate==max(split_times_ES$point_estimate)),]

color_mainSpain<-c("green4", "#6A3D9A", "maroon", "dodgerblue")
tmp2<-coal[coal$group1=="ES03" & coal$group2=="FR01",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df = 20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df = 20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])


rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E3F1_rCCRtab<- rCCRtab
E3F1_tmp2<- tmp2
E3F1_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color=color_mainSpain[1])+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR), color=color_mainSpain[1])+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2, fill= color_mainSpain[1])+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="Relative Cross Coalescent Rate", x= "", title="ES03 vs FR01")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)



tmp2<-coal[coal$group1=="ES04" & coal$group2=="FR01",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)

E4F1_tmp2 <- tmp2
E4F1_rCCRtab <- rCCRtab
E4F1_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point( color=color_mainSpain[2])+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR), color=color_mainSpain[2])+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2, fill=color_mainSpain[2])+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "", title="ES04 vs FR01")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)


tmp2<-coal[coal$group1=="ES17" & coal$group2=="FR01",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)

E17F1_tmp2<- tmp2
E17F1_rCCRtab <- rCCRtab
E17F1_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color=color_mainSpain[3])+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color=color_mainSpain[3])+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill=color_mainSpain[3])+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "", title="ES17 vs FR01")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)

tmp2<-coal[coal$group1=="ES25" & coal$group2=="FR01",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)

E25F1_tmp2<- tmp2
E25F1_rCCRtab <- rCCRtab
E25F1_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color=color_mainSpain[3])+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color=color_mainSpain[3])+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill=color_mainSpain[3])+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "", title="ES25 vs FR01")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)
####
E3F1_rCCRtab$pop <- "ES03"
E4F1_rCCRtab$pop <- "ES04"
E17F1_rCCRtab$pop <- "ES17"
E25F1_rCCRtab$pop <- "ES25"
E3F1_tmp2$pop <- "ES03"
E4F1_tmp2$pop <- "ES04"
E17F1_tmp2$pop <- "ES17"
E25F1_tmp2$pop <- "ES25"


EF1_tmp2<-rbind(E3F1_tmp2, E4F1_tmp2,E17F1_tmp2,E25F1_tmp2 )
EF1_rCCRtab<-rbind(E3F1_rCCRtab, E4F1_rCCRtab,E17F1_rCCRtab, E25F1_rCCRtab)


EF1_rCCR<-ggplot(data=EF1_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EF1_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EF1_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.22, 0.75), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES25", "ES17", "ES04", "ES03"), values = color_mainSpain[4:1], labels=c(
    expression("FR01 vs ES25, " ~ T[div] ~ "= 200 ka"),
    expression("FR01 vs ES17, " ~ T[div] ~ "= 215 ka"),
    expression("FR01 vs ES04, " ~ T[div] ~ "= 229 ka"),
    expression("FR01 vs ES03, " ~ T[div] ~ "= 231 ka")
  ))+
  scale_color_manual(breaks = c("ES25", "ES17", "ES04", "ES03"), values = color_mainSpain[4:1], labels=c(
    expression("FR01 vs ES25, " ~ T[div] ~ "= 200 ka"),
    expression("FR01 vs ES17, " ~ T[div] ~ "= 215 ka"),
    expression("FR01 vs ES04, " ~ T[div] ~ "= 229 ka"),
    expression("FR01 vs ES03, " ~ T[div] ~ "= 231 ka")
  ))
  # annotate(geom = "text", label= "Split time = 197 ka - 245 ka ", x = 5*10^4, y = 1.15, size=10)

EF1_rCCR<-ggplot(data=EF1_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EF1_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EF1_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.22, 0.75), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES25", "ES17", "ES04", "ES03"), values = color_mainSpain[4:1], labels=c(
    expression("FR01 vs ES25"),
    expression("FR01 vs ES17"),
    expression("FR01 vs ES04"),
    expression("FR01 vs ES03")
  ))+
  scale_color_manual(breaks = c("ES25", "ES17", "ES04", "ES03"), values = color_mainSpain[4:1], labels=c(
    expression("FR01 vs ES25"),
    expression("FR01 vs ES17"),
    expression("FR01 vs ES04"),
    expression("FR01 vs ES03")
  ))

##within CAN
tmp2<-coal[coal$group1=="ES03" & coal$group2=="ES04",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E3E4_rCCRtab<- rCCRtab
E3E4_rCCRtab$pop <- "ES03 vs ES04"
E3E4_tmp2 <- tmp2
E3E4_tmp2$pop <- "ES03 vs ES04"

E3E4_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#4F5F8F")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#4F5F8F")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#4F5F8F")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="Relative Cross Coalescent Rate", x= "Years ago", title="ES03 vs ES04")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)


tmp2<-coal[coal$group1=="ES03" & coal$group2=="ES17",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E3E17_rCCRtab <- rCCRtab
E3E17_rCCRtab$pop <- "ES03 vs ES17"
E3E17_tmp2 <- tmp2
E3E17_tmp2$pop <- "ES03 vs ES17"
E3E17_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#804000")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#804000")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#804000")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "Years ago", title="ES03 vs ES17")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)


tmp2<-coal[coal$group1=="ES03" & coal$group2=="ES25",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E3E25_rCCRtab <- rCCRtab
E3E25_rCCRtab$pop <- "ES03 vs ES25"
E3E25_tmp2 <- tmp2
E3E25_tmp2$pop <- "ES03 vs ES25"
E3E25_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#0F8E80")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#0F8E80")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#0F8E80")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "Years ago", title="ES03 vs ES17")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)



tmp2<-coal[coal$group1=="ES04" & coal$group2=="ES17",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E4E17_rCCRtab <- rCCRtab
E4E17_rCCRtab$pop <- "ES04 vs ES17"
E4E17_tmp2 <- tmp2
E4E17_tmp2$pop <- "ES04 vs ES17"

E4E17_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#85306D")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#85306D")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#85306D")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "Years ago", title="ES04 vs ES17")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)


tmp2<-coal[coal$group1=="ES04" & coal$group2=="ES25",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E4E25_rCCRtab <- rCCRtab
E4E25_rCCRtab$pop <- "ES04 vs ES25"
E4E25_tmp2 <- tmp2
E4E25_tmp2$pop <- "ES04 vs ES25"
E4E25_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#4466CC")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#4466CC")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#4466CC")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "Years ago", title="ES04 vs ES25")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)
color_mainSpain


tmp2<-coal[coal$group1=="ES17" & coal$group2=="ES25",]
tmp2<- tmp2[!is.na(tmp2$rCCR) &tmp2$years!=0 ,]
log_years<- log(tmp2$years)
spline_fit <- smooth.spline(log_years, tmp2$rCCR, df = 20)  # Adjust df for more or less flexibility
spline_fit_min <- smooth.spline(log_years, tmp2$rCCR.min, df=20)  # Adjust df for more or less flexibility
spline_fit_max <- smooth.spline(log_years, tmp2$rCCR.max, df=20)  # Adjust df for more or less flexibility
max_years<- max(log_years)
min_years<- min(log_years)

new_years<- seq(min_years,max_years, length.out=1000)
predicted_rCCR <- predict(spline_fit,x=new_years)
predicted_rCCR_min <- predict(spline_fit_min, x=new_years)
predicted_rCCR_max <- predict(spline_fit_max, x=new_years)

plot(log_years, tmp2$rCCR, pch = 16, main = "Spline Fit", xlab = "Years", ylab = "rCCR")
lines(new_years, predicted_rCCR$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_min$y, col = "green", lwd = 2)
lines(new_years, predicted_rCCR_max$y, col = "green", lwd = 2)

closest_index <- min(which(abs(predicted_rCCR$y - 0.5)<0.01))
closest_index_min <- min(which(abs(predicted_rCCR_min$y - 0.5)<0.01))
closest_index_max <- min(which(abs(predicted_rCCR_max$y - 0.5)<0.01))

exp(predicted_rCCR$x[closest_index])
exp(predicted_rCCR$x[closest_index_min])
exp(predicted_rCCR$x[closest_index_max])  

rCCRtab<-data.frame(years=new_years, rCCR=predicted_rCCR$y, rCCRmin=predicted_rCCR_max$y, rCCRmax=predicted_rCCR_min$y)
E17E25_rCCRtab <- rCCRtab
E17E25_rCCRtab$pop <- "ES17 vs ES25"
E17E25_tmp2 <- tmp2
E17E25_tmp2$pop <- "ES17 vs ES25"
E17E25_rCCR<-ggplot(data=tmp2, aes(x = years,y= rCCR))+
  geom_point(color="#4F4880")+
  geom_line(data = rCCRtab, aes(x=exp(years), y=rCCR),color="#4F4880")+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax), alpha=0.2,fill="#4F4880")+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.7, 0.7), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="", x= "Years ago", title="ES04 vs ES25")+
  annotate(geom = "text", label=paste0("Split time = ", as.character(round(exp(predicted_rCCR$x[closest_index])/1000))," ka (", as.character(round(exp(predicted_rCCR$x[closest_index_max])/1000)), " ka - ",as.character(round(exp(predicted_rCCR$x[closest_index_min])/1000))," ka)"), x = 5*10^4, y = 1.15, size=10)


EE_tmp2<-rbind(E3E4_tmp2, E3E17_tmp2,E4E17_tmp2,E3E25_tmp2, E4E25_tmp2, E17E25_tmp2 )
EE_rCCRtab<-rbind(E3E4_rCCRtab, E3E17_rCCRtab,E4E17_rCCRtab, E3E25_rCCRtab, E4E25_rCCRtab,E17E25_rCCRtab)

EE_rCCR<-ggplot(data=EE_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EE_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EE_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.22, 0.75), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
    values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
    expression("ES03 vs ES04, " ~ T[div] ~ "= 19 ka"),
    expression("ES04 vs ES17, " ~ T[div] ~ "= 70 ka"),
    expression("ES03 vs ES17, " ~ T[div] ~ "= 78 ka"),
    expression("ES17 vs ES25, " ~ T[div] ~ "= 121 ka"),
    expression("ES04 vs ES25, " ~ T[div] ~ "= 140 ka"),
    expression("ES03 vs ES25, " ~ T[div] ~ "= 143 ka")
  ))+
  scale_color_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
    values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
    expression("ES03 vs ES04, " ~ T[div] ~ "= 19 ka"),
    expression("ES04 vs ES17, " ~ T[div] ~ "= 70 ka"),
    expression("ES03 vs ES17, " ~ T[div] ~ "= 78 ka"),
    expression("ES17 vs ES25, " ~ T[div] ~ "= 121 ka"),
    expression("ES04 vs ES25, " ~ T[div] ~ "= 140 ka"),
    expression("ES03 vs ES25, " ~ T[div] ~ "= 143 ka")
  ))


EE_rCCR<-ggplot(data=EE_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EE_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EE_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = 0, axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                    values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                      expression("ES03 vs ES04, " ~ T[div] ~ "= 19 ka"),
                      expression("ES04 vs ES17, " ~ T[div] ~ "= 70 ka"),
                      expression("ES03 vs ES17, " ~ T[div] ~ "= 78 ka"),
                      expression("ES17 vs ES25, " ~ T[div] ~ "= 121 ka"),
                      expression("ES04 vs ES25, " ~ T[div] ~ "= 140 ka"),
                      expression("ES03 vs ES25, " ~ T[div] ~ "= 143 ka")
                    ))+
  scale_color_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                     values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                       expression("ES03 vs ES04, " ~ T[div] ~ "= 19 ka"),
                       expression("ES04 vs ES17, " ~ T[div] ~ "= 70 ka"),
                       expression("ES03 vs ES17, " ~ T[div] ~ "= 78 ka"),
                       expression("ES17 vs ES25, " ~ T[div] ~ "= 121 ka"),
                       expression("ES04 vs ES25, " ~ T[div] ~ "= 140 ka"),
                       expression("ES03 vs ES25, " ~ T[div] ~ "= 143 ka")
                     ))

EE_rCCR<-ggplot(data=EE_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EE_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EE_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = c(0.22, 0.75), axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                    values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                      expression("ES03 vs ES04"),
                      expression("ES04 vs ES17"),
                      expression("ES03 vs ES17"),
                      expression("ES17 vs ES25"),
                      expression("ES04 vs ES25"),
                      expression("ES03 vs ES25")
                    ))+
  scale_color_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                     values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                       expression("ES03 vs ES04"),
                       expression("ES04 vs ES17"),
                       expression("ES03 vs ES17"),
                       expression("ES17 vs ES25"),
                       expression("ES04 vs ES25"),
                       expression("ES03 vs ES25")
                     ))

EE_rCCR<-ggplot(data=EE_tmp2, aes(x = years,y= rCCR))+
  geom_point(aes(color=pop), size=3)+
  geom_line(data = EE_rCCRtab, aes(x=exp(years), y=rCCR, color=pop), size=1.2)+
  scale_y_continuous(limits = c(0,1.2), breaks=seq(0,1,by=0.5))+
  geom_ribbon(data = EE_rCCRtab, aes(x=exp(years), ymin=rCCRmin, ymax=rCCRmax, fill=pop), alpha=0.2)+
  scale_x_continuous(trans = log_trans(10), limits = c(1e3,2682690), breaks = c(1e3, 1e4,1e5,1e6), labels = c(expression(paste(10^3)),expression(paste(10^4)),expression(paste(10^5)),expression(paste(10^6))), expand = c(0,0)) +
  annotation_logticks(sides = "tb", long = unit(0.5, "cm"), mid = unit(0.25, "cm"), short = unit(0.125, "cm")) +
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, colour = "white"), panel.border = element_rect(fill=NA, color="black"), legend.position = "none", axis.title=element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 22), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5))+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  labs(y="rCCR", x= "Years ago", title="", color="", fill="")+
  scale_fill_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                    values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                      expression("ES03 vs ES04"),
                      expression("ES04 vs ES17"),
                      expression("ES03 vs ES17"),
                      expression("ES17 vs ES25"),
                      expression("ES04 vs ES25"),
                      expression("ES03 vs ES25")
                    ))+
  scale_color_manual(breaks = c("ES03 vs ES04", "ES04 vs ES17", "ES03 vs ES17", "ES17 vs ES25", "ES04 vs ES25", "ES03 vs ES25"),
                     values = c("#4F5F8F","#85306D","#804000", "#4F4880","#4466CC" , "#0F8E80"), labels=c(
                       expression("ES03 vs ES04"),
                       expression("ES04 vs ES17"),
                       expression("ES03 vs ES17"),
                       expression("ES17 vs ES25"),
                       expression("ES04 vs ES25"),
                       expression("ES03 vs ES25")
                     ))
library(cowplot)


rCCR_plot<-plot_grid(plotlist = list(EF1_rCCR, EE_rCCR ), cols = 2)
plot_grid(plotlist = list(popsize, rCCR_plot ), cols = 1)
plot_grid(plotlist = list(E3F1_rCCR,E4F1_rCCR,E17F1_rCCR, E3E4_rCCR, E3E17_rCCR, E4E17_rCCR), cols = 3)




##################correlation analysis split times

##initial split time anylsis
split_times_ES<- split_times[grep(x = split_times$pop1,pattern = "ES")[grep(x = split_times$pop1,pattern = "ES") %in% grep(x = split_times$pop2,pattern = "ES")],]
# split_times_ES<- split_times_ES[order(split_times_ES$point_estimate),]
i<-1
while(i <nrow(split_times_ES)){
  # split_times_ES<- 
  split_times_ES<-split_times_ES[!(split_times_ES$pop1==split_times_ES$pop2[i] & split_times_ES$pop2==split_times_ES$pop1[i]),]
  i<-i+1
}
split_times_ES


coords<-read.table("~/Data/Spanish_adaptation/coords.txt", header=T, sep="\t")
coords<- coords[,1:3]
split_times_ES$max_lon <- NA
for(i in 1:nrow(split_times_ES)){
  split_times_ES$max_lon[i] <- max(coords$lon[split_times_ES$pop1[i]==coords$group],coords$lon[split_times_ES$pop2[i]==coords$group])  
}
split_times_ES$point_estimate <- as.numeric(split_times_ES$point_estimate)
cor(x = split_times_ES$max_lon,y = split_times_ES$point_estimate)
cor.test(x = split_times_ES$max_lon,y = split_times_ES$point_estimate)


ggplot(data = split_times_ES, aes(x=max_lon, y = point_estimate))+
  geom_point(color="black")+
  labs(x="Longitude", y=expression(T[div] ~ "[x" * 10^3 *"]"))+
  scale_x_continuous(limits=c(-7,-3), expand = c(0.1,0))+
  scale_y_continuous(limits=c(0, 150000), breaks = seq(0,150000,by=50000), labels = seq(0,150,by=50), expand = c(0.1,0))+
  theme(panel.background = element_blank(),axis.text = element_text(size=40), axis.title = element_text(size=40), axis.ticks.length = unit(1.5, "cm"))+
  geom_segment(x = -7, xend=-3, y=-15000, yend=-15000)+
  geom_segment(x = -7.4, xend=-7.4, y=0, yend=150000)

ggplot(data = split_times_ES, aes(x = max_lon, y = point_estimate)) +
  geom_point(color = "black", size=3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add linear trend line
  labs(
    x = "Longitude", 
    # y = expression(T[div] ~ "[" * 10^3 * "]")
    # y = expression("Years ago [" * 10^3 * "]" ~ "(" * T[div] * ")")
    y = expression("Years ago [" * 10^3 * "]")
  ) +
  scale_x_continuous(
    limits = c(-7, -3), 
    breaks = seq(-7, -3, by = 1),  # Set breaks at each integer longitude
    labels = function(x) paste0(abs(x), "° W"),  # Convert to degrees west format
    expand = c(0.1, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 150000), 
    breaks = seq(0, 150000, by = 50000), 
    labels = seq(0, 150, by = 50), 
    expand = c(0.1, 0)
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 40), 
    axis.title = element_text(size = 40), 
    axis.ticks.length = unit(1.5, "cm")
  ) +
  geom_segment(x = -7, xend = -3, y = -15000, yend = -15000) +
  geom_segment(x = -7.4, xend = -7.4, y = 0, yend = 150000) +
  annotate(
    "text", 
    x = -6,  # Middle of the x-axis
    y = 125000,  # Middle of the y-axis
    label = "r = 0.836, p-value < 2.2e-16", 
    size = 10
  )

  

########after revision with permutations###
library(geosphere)

split_times_ES<- split_times[grep(x = split_times$pop1,pattern = "ES")[grep(x = split_times$pop1,pattern = "ES") %in% grep(x = split_times$pop2,pattern = "ES")],]
i<-1
while(i <nrow(split_times_ES)){
  split_times_ES<-split_times_ES[!(split_times_ES$pop1==split_times_ES$pop2[i] & split_times_ES$pop2==split_times_ES$pop1[i]),]
  i<-i+1
}
split_times_ES


coords<-read.table("~/Data/Spanish_adaptation/coords.txt", header=T, sep="\t")
coords<- coords[,1:3]
split_times_ES$lon1 <-NA
split_times_ES$lon2 <-NA
split_times_ES$lat1 <-NA
split_times_ES$lat2 <-NA

split_times_ES$max_lon <- NA
split_times_ES$min_lon <- NA
split_times_ES$mid_lon <- NA
split_times_ES$AvLonDiff_ToMaxDiff<- NA
for(i in 1:nrow(split_times_ES)){

  lon1 <- coords$lon[split_times_ES$pop1[i] == coords$group]
  lon2 <- coords$lon[split_times_ES$pop2[i] == coords$group]
  lat1 <- coords$lat[split_times_ES$pop1[i] == coords$group]
  lat2 <- coords$lat[split_times_ES$pop2[i] == coords$group]

  split_times_ES$lon1[i] <-lon1
  split_times_ES$lon2[i] <-lon2
  split_times_ES$lat1[i] <-lat1
  split_times_ES$lat2[i] <-lat2

  # reorder so pop1 = western, pop2 = eastern
  if(lon1 > lon2){

    tmp <- split_times_ES$pop1[i]
    split_times_ES$pop1[i] <- split_times_ES$pop2[i]
    split_times_ES$pop2[i] <- tmp

    tmp_lon <- lon1
    lon1 <- lon2
    lon2 <- tmp_lon
  }

  split_times_ES$max_lon[i] <- lon2
  split_times_ES$min_lon[i] <- lon1
  split_times_ES$mid_lon[i] <- (lon1+lon2)/2
}
split_times_ES$point_estimate <- as.numeric(split_times_ES$point_estimate)

max_lon <- max(split_times_ES$lon1, split_times_ES$lon2)
split_times_ES$AvLonDiff_ToMaxDiff <-  (max_lon - split_times_ES$lon1)/2 + (max_lon - split_times_ES$lon2)/2
split_times_ES$MinLonDiff_ToMaxDiff <-  pmin((max_lon - split_times_ES$lon1),(max_lon - split_times_ES$lon2))




# observed correlation
#geographic distance
obs_corr_dist <- cor(split_times_ES$dist,
               split_times_ES$point_estimate)
#mid point longitude of each pair
obs_corr_mid <- cor(split_times_ES$mid_lon,
                split_times_ES$point_estimate)
#max (most eastern) longitude of each pair
obs_corr_max <- cor(split_times_ES$max_lon,
                split_times_ES$point_estimate)
#min (most western) longitude of each pair
obs_corr_min <- cor(split_times_ES$min_lon,
                split_times_ES$point_estimate)
#minimal geographic distance to the most eastern population in the set
obs_corr_Mindist2maxlonpop <- cor(split_times_ES$Mindist2maxlonpop,
                                  split_times_ES$point_estimate)
#average geographic distance of the population pair to the most eastern population in the set
obs_corr_Meandist2maxlonpop <-cor(split_times_ES$Meandist2maxlonpop
                                  , split_times_ES$point_estimate)


####permutations randomly shuffle coordinates between popultions

set.seed(1)

nperm <- 1000000

perm_corr_mid <- numeric(nperm)
perm_corr_max <- numeric(nperm)
perm_corr_min <- numeric(nperm)
perm_corr_dist <- numeric(nperm)
perm_corr_Mindist2maxlonpop <- numeric(nperm)
perm_corr_Meandist2maxlonpop <- numeric(nperm)


length(unique(split_times_ES$pop1))
length(pop_lons)

pop_coord <- coords[,c(2,3)]
row.names(pop_coord) <- coords$group
pop_coord<- pop_coord[rownames(pop_coord)!="ES15",]



pop_names <- rownames(pop_coord)

i1 <- match(split_times_ES$pop1, pop_names)
i2 <- match(split_times_ES$pop2, pop_names)

lon <- pop_coord[,1]
lat <- pop_coord[,2]
y <- split_times_ES$point_estimate

dist_mat <- distm(
  x = cbind(lon, lat),
  fun = distHaversine
)

dist2maxlonpop <- dist_mat[13,]

op<-0.01*nperm
i<-1
for(i in 1:nperm){

  idx <- sample(seq_along(lon))

  plon <- lon[idx]

  perm_mid_lon <- (plon[i1] + plon[i2]) / 2
  perm_max_lon <- pmax(plon[i1], plon[i2])
  perm_min_lon <- pmin(plon[i1], plon[i2])

  perm_dist <- dist_mat[cbind(idx[i1], idx[i2])]
  perm_dist2maxlonpop <- dist2maxlonpop[idx]
  perm_Mindist2maxlon <- pmin(perm_dist2maxlonpop[i1],perm_dist2maxlonpop[i2])
  perm_Meandist2maxlon <- (perm_dist2maxlonpop[i1])/2 + (perm_dist2maxlonpop[i2])/2

  perm_corr_dist[i] <- cor(perm_dist, y)
  perm_corr_mid[i]  <- cor(perm_mid_lon, y)
  perm_corr_max[i]  <- cor(perm_max_lon, y)
  perm_corr_min[i]  <- cor(perm_min_lon, y)
  perm_corr_Mindist2maxlonpop[i] <- cor(perm_Mindist2maxlon, y)
  perm_corr_Meandist2maxlonpop[i] <- cor(perm_Meandist2maxlon, y)

  if(i%%op==0){
    print(floor(100*i/nperm))
  }
}

p_perm_dist <- mean(abs(perm_corr_dist) >= abs(obs_corr_dist))
p_perm_mid <- mean(abs(perm_corr_mid) >= abs(obs_corr_mid))
p_perm_max <- mean(abs(perm_corr_max) >= abs(obs_corr_max))
p_perm_min <- mean(abs(perm_corr_min) >= abs(obs_corr_min))
p_perm_Mindist2maxlonpop <- mean(abs(perm_corr_Mindist2maxlonpop) >= abs(obs_corr_Mindist2maxlonpop))
p_perm_Meanndist2maxlonpop <- mean(abs(perm_corr_Meandist2maxlonpop) >= abs(obs_corr_Meandist2maxlonpop))

obs_corr_dist
#0.8535486
p_perm_dist
#1e-05

obs_corr_mid
#0.6866079
p_perm_mid
#0.000256

obs_corr_max
#0.8360898
p_perm_max
#2.7e-05


obs_corr_min
#0.1404573
p_perm_min
#0.470391

obs_corr_Mindist2maxlonpop
#-0.8371181
p_perm_Mindist2maxlonpop
#2.6e-05

obs_corr_Meandist2maxlonpop
#-0.6880119
p_perm_Meanndist2maxlonpop
#0.000238
