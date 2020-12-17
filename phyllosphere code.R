#phyllosphere paper 

leaf.asv<- read_biom('no-mito-chloro.biom')
leaf.tax<-read.csv('no.chloro.mito.leaf.tax.separated.csv')
leaf.sam<-read.csv('leaf.sample.id.csv')

lew<- phyloseq(otu_table(leaf.asv, taxa_are_rows = TRUE), sample_data(leaf.sam), tax_table(leaf.tax)) 
lew<-subset_samples(leaf.otu, site== "Llewelyn")
lew.1<- prune_samples(sample_sums(lew) > 0, lew)

#diversity 
inv.simp<-estimate_richness(leaf.otu, measures='InvSimpson') #indication of the richness in a community with uniform evenness, gives weight to common species 

shannon<- estimate_richness(leaf.otu, measures='Shannon') #combines evenness and richness, efected by sampling effort 

ace<-estimate_richness(leaf.otu, measures='ACE') #uses # of species, sorts into abundant vs not, then only does prescense/absence 

#species richness
abun <- t(data.frame(otu_table(lew.1)))

spec<-as.data.frame(specnumber(abun))

ggplot(lew.mean, aes(x=week, y=ASVs, group= 1 )) + geom_point(size=4) + geom_line() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), panel.background = element_rect(fill='white')) + labs(x='Week', y='# ASVs') +
  ggtitle("Llewellyn Observed Diversity by Week") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.x = element_text(size=20), axis.text.y  = element_text(size=12))


#PERMANOVAS
L.weight<-phyloseq::distance(lew.1, method='wunifrac')
adonis(L.weight~stage*NAM, data= L.sam) 


L.unweight<-phyloseq::distance(lew.1, method='unifrac')
adonis(L.unweight~stage*NAM, data=L.sam) 


#winnowed community = Expanded Core
#win.weight, win.unweight, win.sam

win.sam<-data.frame(sample_data(L.win7))
win.wt<-as.matrix(win.weight)
win.unwt<-as.matrix(win.unweight)

rda5<-capscale(win.wt~NAM+stage, data=win.sam)
rda5
plot.new()
plot(rda5, display=c('cn'), scaling='species')
points(rda5, pch=20, col = L.sam$stage)
text(rda5, dis="cn", lwd = 2, cex = 0.00000001)

anova(rda5, by='terms') 


rda6<-capscale(win.unwt~NAM*stage, data=win.sam)
rda6

plot.new()
plot(rda6, display=c('bp'), scaling='species')
points(rda6, pch=20, col = L.sam$stage)
text(rda6, dis="bp", lwd = 2, cex = 0.00000001)

anova(rda6, by='terms') 

#permanovas, not constrained
per1<-adonis(L.wt1~stage*NAM*block, data=l.sam)
per1

per2<-adonis(L.unwt1~stage*NAM*block, data=l.sam)
per2

per3<-adonis(win.wt~stage*NAM*block, data=win.sam.1)
per3

per4<-adonis(win.unwt~stage*NAM*block, data=win.sam.1)
per4


### just do rdas w/ stage, then extract residulas then do rda on those w/ NAM

rda7<-capscale(L.wt1~stage, data=L.sam)
rda7
res.1<-as.data.frame(rstandard(rda7, type= "canoco"))

c.res1<-capscale(res.1~NAM, data=L.sam)
plot(c.res1)
per5<-adonis(res.1~NAM, data=L.sam)
per5

c.res1
anova(c.res1)


plot(c.res1, display=c('bp', 'wa'), scaling='species')
points(c.res1, pch=20, col = L.sam$NAM)
text(rda7, dis="bp", lwd = 2, cex = 0.00000001)


#whole unweighted
rda8<-capscale(L.unwt1~stage, data=L.sam)
rda8
res.2<-as.data.frame(rstandard(rda8, type= "canoco"))

c.res2<-capscale(res.2~NAM, data=L.sam)
per6<-adonis(res.2~NAM, data=L.sam)
per6

anova(c.res2, by='terms')

plot(c.res2, display=c('bp', 'wa'), scaling='species')
points(c.res2, pch=20, col = L.sam$NAM)

#win weighted
rda9<-capscale(win.wt~stage, data=win.sam)
rda9
res.3<-as.data.frame(rstandard(rda9, type= "canoco"))

c.res3<-capscale(res.3~NAM, data=win.sam)
per7<-adonis(res.3~NAM, data=win.sam)
per7

anova(c.res3)

plot(c.res3, display=c('bp', 'wa'), scaling='species')
points(c.res3, pch=20, col = win.sam$NAM)

#win unweighted
rda10<-capscale(win.unwt~stage, data=win.sam)
rda10
res.4<-as.data.frame(rstandard(rda10, type= "canoco"))

c.res4<-capscale(res.4~NAM, data=win.sam)
per8<-adonis(res.4~NAM, data=win.sam)
per8

c.res4

anova(c.res4,by='terms')

plot(c.res4, display=c('bp', 'wa'), scaling='species')
points(c.res4, pch=20, col = win.sam$NAM)


#do dbRDA of Flower and Pod

pre.2<-subset_samples(lew.1, week == c("1", '2', '3', '4'))
post.5<-subset_samples(lew.1, week == c("5", '6', '7', '8', '9', '10'))

pre.uni.dist<-phyloseq::distance(pre.2,   "wunifrac")
unpre.uni.dist<-phyloseq::distance(pre.2,   "uunifrac")
pre.samp<-data.frame(sample_data(pre.2))


post.uni.dist<-phyloseq::distance(post.5,   "wunifrac")
unpost.uni.dist<-phyloseq::distance(post.5,   "uunifrac")
#already have post.samp

pre.wt.rda<-capscale(pre.uni.dist~NAM*stage, data=pre.samp)
anova(pre.wt.rda, by='terms')

per9<-adonis(pre.uni.dist~stage*NAM*block, data=pre.sam)
per9


pre.unwt.rda<-capscale(unpre.uni.dist~NAM*stage, data=pre.samp)
pre.unwt.rda
per10<-adonis(unpre.uni.dist~stage*NAM*block, data=pre.sam)
per10


post.wt.rda<-capscale(post.uni.dist~NAM*stage, data=post.samp)
post.wt.rda
per11<-adonis(post.uni.dist~stage*NAM*block, data=post.sam)
per11


post.unwt.rda<-capscale(unpost.uni.dist~NAM*stage, data=post.samp)
post.unwt.rda
per12<-adonis(unpost.uni.dist~stage*NAM*block, data=post.sam)
per12
