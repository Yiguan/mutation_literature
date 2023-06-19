
library('readxl')
library('data.table')
library('dplyr')
library('tidyr')

library('ggplot2')
library("cowplot")
library("RColorBrewer")
library("grid")
library("gtable")
library('ggpubr')
library("ggrepel")
library("ggtree")

library("MCMCglmm")
library("Hmisc")


mycol <- brewer.pal(12, 'Set3')

#human_col       = mycol[4]
primates_col    = mycol[4]
mammals_col     = mycol[6]
birds_col       = mycol[1]
reptiles_col   = mycol[5]
fish_col        = mycol[7]
arthropods_col  = mycol[11]
nematodes_col   = mycol[12]
plants_col      = mycol[10]
fungi_col      = mycol[8]
others_col      = mycol[9]


###############################################
## FUNCTIONS

# infer effective number of events based on interval assuming binomial distribution
inferEvents <- function(um, u_lower, u_upper){
    delta = 100
    i = 1
    while (TRUE){
        tt <- binconf(i,i/um)
        t_lower <- tt[2]
        t_upper <- tt[3]
        low_diff <- (t_lower - u_lower)/u_lower
        up_diff <- (t_upper - u_upper)/u_upper
        tmp <- abs(low_diff) + abs(up_diff)
        if(tmp > delta){
            break
        }else{
            delta <- tmp
            i = i + 1
        }
    }
    return(i-1)
}

# get lower CI
inferLowerCI <- function(x,n){
    tt <- binconf(x,n)
    return(tt[2])
}

# get upper CI
inferUpperCI <- function(x,n){
    tt <- binconf(x,n)
    return(tt[3])   
}

################################################################################

bb <- read_excel('mutation_rate_literature_updating3.xlsx', sheet='raw') %>% 
    filter(Method %in% c('MA','PO')) %>%
    select(-c('range_low','range_up','SampleSize','Generations','Sequencing',
              'Tools', 'Depth', 'Citation2', 'Link', 'Populations')) %>%
    as.data.table() %>%
    filter(! PID %in% c(2023002,2023003,2023004))
#str(bb)

bb[,PID] %>% unique() %>% length() # Number of studies
bb$Species %>% unique() %>% length() # Mumber of species

###############################################################################

b_year <- bb %>% filter(Date<2023) %>% 
    select(c('PID', 'Method')) %>% unique() %>%
    tidyr::extract(PID, 'Year', '(....).*') %>%
    group_by(Year, Method) %>% count() %>%
    mutate(Year = as.integer(Year))

b_year$Method <- factor(b_year$Method, levels=c('MA', 'PO'))


# bar plot for number of studies by year
P_bar <- ggplot(b_year, aes(Year, n, fill=Method)) +
    geom_bar(stat='identity') +
    scale_y_continuous(name='Number of Studies',expand = c(0,0), limits = c(0,20)) +
    scale_x_continuous(breaks = seq(min(b_year$Year), max(b_year$Year), by = 1)) +
    scale_fill_manual(name='', values = c('MA'='#009999','PO'='darkorange3')) +
    theme_classic() +
    guides(fill=guide_legend(nrow=2)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4),
          #panel.grid.major = element_line(),
          axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          legend.text = element_text(size=16),
          legend.position = c(0.2,0.8),
          legend.key.size = unit(2, 'lines'),
          plot.margin = margin(3, 0.4, 3, 0.4, "cm"))


################################################################################

b_species <- bb %>% filter(Date<2023) %>% select(c('PID', 'Group2')) %>% 
    unique() %>%
    group_by(Group2) %>% count()

b_species$Group2 <- factor(b_species$Group2, 
                          levels=c('primates','mammals','birds','reptiles',
                                   'fish','arthropods','nematodes','plants',
                                   'fungi','others'))

b_species <- b_species %>% arrange(desc(Group2)) %>%
    mutate(prop = n/sum(b_species$n) * 100) %>% as.data.frame() %>%
    mutate(ypos = cumsum(prop) - 0.5*prop) %>%
    mutate(ll = sprintf('%s\n%.1f%%', Group2, prop))


scol <- c('primates'=primates_col,'mammals'=mammals_col, 'birds'=birds_col, 
          'reptiles'=reptiles_col, 'fish'=fish_col,'arthropods'=arthropods_col,
          'nematodes'=nematodes_col,'plants'=plants_col,
          'fungi'=fungi_col,'others'=others_col)


# Pie chart for proportion of species
P_pie_species <- ggplot(b_species, aes(x='', y=prop, fill=Group2)) +
    geom_bar(width = 1, stat='identity', color='white') +
    scale_fill_manual(values=scol) +
    coord_polar('y',start=0) +
    theme_void() +
    #geom_text(aes(y=ypos, label=Group), size=5)+ 
    theme(legend.position = 'none') +
    geom_label_repel(aes(y=ypos, label=ll), size=4.5, nudge_x = 0.3, 
                     min.segment.length=unit(6,'lines'),
                     color='black', box.padding = 0.5, fill='white')

################################################################################

p12 <- plot_grid(P_pie_species, P_bar, nrow=1, labels = c( '(A)', '(B)'), label_x = 0.1)

# ggsave('study_summary.jpg', p12, width=12, height=8, units = 'in', dpi = 300)
ggsave('study_summary.pdf', p12, width=12, height=8, units = 'in', dpi = 300)

##
# b1 <- bb %>% select(c('PID', 'Group', 'Method')) %>% unique()
# b1[Group!='human', Method] %>% table
# b_year %>% as.data.table %>% filter(Year>=2010) %>%
#     select(n) %>% sum()


################################################################################
###### PART TWO
################################################################################



bb_snp <- bb %>% extract(CITATION, 'Author', '^(.*?),') %>%
    tidyr::unite('PID2','Author','Date', sep=',') %>%
    filter(TYPE=='snp') %>% 
    filter(Unit=='bp/generation') %>%
    mutate(u_mean = as.numeric(u_mean)) %>% 
    mutate(u_lower = as.numeric(u_lower)) %>%
    mutate(u_upper = as.numeric(u_upper)) %>%
    mutate(SE = as.numeric(SE)) %>%
    mutate(Events = as.integer(Events)) %>%
    filter(!Comments %like% 'exome') %>%
    filter(!Comments %like% 'male') %>% 
    filter(!Comments %like% 'sperm') %>% 
    select(c('PID','PID2','SID','u_mean','u_lower','u_upper','Events','Callable', 'SE',
             'reproduction','Species','Method', 'Group2')) %>%
    as.data.table

# Group close strain togther for phylogeny

bb_snp[Species %in% c('Anopheles coluzzii', 'Anopheles stephensi'), Species:='Anopheles gambiae']
bb_snp[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species:='Ostreococcus tauri']
bb_snp[Species == 'Chlamydomonas incerta', Species:='Chlamydomonas reinhardtii']
bb_snp[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), Species:='Paramecium tetraurelia']
bb_snp[Species =='Saccharomycodes ludwigii', Species:='Hanseniaspora valbyensis']

# Calculate interval if SE and u_mean are provided, u +/- 1.96 SE ?

bb_snp[is.na(Events) & is.na(u_lower) & !is.na(SE) & !is.na(u_mean), 
       u_lower:=u_mean-1.96*SE]
bb_snp[is.na(Events) & is.na(u_upper) & !is.na(SE) & !is.na(u_mean), 
       u_upper:=u_mean+1.96*SE]

# Estimate Events based on mean and interval
bb_snp[is.na(Events) & !is.na(u_lower), 
       Events:=apply(.SD,1,function(x) inferEvents(x[1],x[2],x[3])), 
       .SDcols=c('u_mean','u_lower','u_upper')]

# Assume 1 Event if interval no available (only u_mean provided)
bb_snp[is.na(Events) & is.na(u_lower), Events:=1]  


bb_snp[, Callable:=Events/u_mean]
bb_snp$PID <- as.character(bb_snp$PID)
############################


bb_snp[,Species2:=Species]
bb_snp[Species=='Amphilophus', Species2:='Amphilophus labiatus']
# bb_snp[Species=='Anopheles stephensi', Species2:='Anopheles gambiae']
# bb_snp[Species=='Anopheles coluzzii', Species2:='Anopheles gambiae']
# bb_snp[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), 
#     Species2:='Paramecium tetraurelia']
#bb_snp[Species=='Chlamydomonas incerta', Species2:='Chlamydomonas reinhardtii']
bb_snp[Species=='Daphnia galeata', Species2:='Daphnia dubia']
#bb_snp[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), 
#       Species2:='Ostreococcus tauri']
#bb_snp[Species=='Saccharomycodes ludwigii', Species2:='Hanseniaspora valbyensis']
# following automatically replaced with neighbouring species by treetime.org
#bb_snp[Species=='Micromonas pusilla', Species2:='Micromonas']
bb_snp[Species=='Picochlorum costavermella', Species2:='Coccomyxa subellipsoidea']
bb_snp[Species=='Rhodotorula toruloides', Species2:='Microbotryomycetes']
bb_snp[Species=='Bombus terrestris', Species2:='Bombus']
bb_snp[Species=='Chironomus riparius', Species2:='Chironomus pallidivittatus']
bb_snp[Species=='Pristionchus pacificus', Species2:='Pristionchus']

bb_snp[Species=='Tupaia chinensis belangeri', Species2:='Tupaia belangeri']
bb_snp[Species=='Ailurus fulgens', Species2:='Ailurus']
bb_snp[Species=='Ceratotherium simum simum', Species2:='Ceratotherium simum']
bb_snp[Species=='Cervus elaphus yarkandensis', Species2:='Cervus hanglu']
bb_snp[Species=='Giraffa camelopardalis', Species2:='Giraffa reticulata']
bb_snp[Species=='Sphaerodactylus inigoi', Species2:='Sphaerodactylus copei']
bb_snp[Species=='Pelecanus crispus', Species2:='Pelecanus occidentalis']
bb_snp[Species=='Emiliania huxleyi', Species2:='Emiliania']
bb_snp[Species=='Phoenicopterus roseus', Species2:='Phoenicopterus']
bb_snp[Species=='Gyps fulvus', Species2:='Gyps']
bb_snp[Species=='Heliconius melpomene', Species2:='Heliconius ethilla']
bb_snp[Species=='Brassica rapa', Species2:='Brassica juncea']
bb_snp[Species=='Oryza sativa', Species2:='Oryza longistaminata']
bb_snp[Species=='Zea mays', Species2:='Zea diploperennis']



bb_snp[Species=='Panthera tigris', Species2:='Panthera tigris tigris']
bb_snp[Species=='Canis lupus', Species2:='Canis lupus lupus']
bb_snp[Species=='Cervus nippon', Species2:='Cervus nippon centralis']
bb_snp[Species=='Cavia aperea', Species2:='Cavia aperea guianae']
bb_snp[Species=='Mus musculus', Species2:='Mus musculus musculus']
bb_snp[Species=='Pan troglodytes', Species2:='Pan troglodytes troglodytes']
bb_snp[Species=='Hylobates lar', Species2:='Hylobates lar lar']
bb_snp[Species=='Macaca mulatta', Species2:='Macaca mulatta vestita']
bb_snp[Species=='Cyanistes caeruleus', Species2:='Cyanistes caeruleus caeruleus']
bb_snp[Species=='Apis mellifera', Species2:='Apis mellifera mellifera']
bb_snp[Species=='Homo sapiens', Species2:='Homo sapiens neanderthalensis']
bb_snp[Species=='Turdus merula', Species2:='Turdus merula merula']
bb_snp[Species=='Gallus gallus domesticus', Species2:='Gallus gallus']

bb_unique_species <- bb_snp[,"Species2"] %>% unique()
fwrite(bb_unique_species, 'filtered_unique_species.txt', col.names = FALSE)

################################################################################


## in bash
## grep ':0.00000000' filtered_unique_species.nwk
## sed -i 's/:0.00000000/:0.01000000/g' filtered_unique_species.nwk




# create phylogeny based on species common names
my_phylo <- read.tree('filtered_unique_species.nwk')
my_phylo <- phytools::force.ultrametric(my_phylo, method='extend')

# rename node labels as '^n\\d+'
my_phylo$node.label <- 
    stringr::str_replace_all(my_phylo$node.label, "'", "") %>% 
    paste0('n', .)

my_phylo <- ladderize(my_phylo)


# rename tip labels as species names
for(nn in my_phylo$tip.label){
    nn1 <- stringr::str_replace_all(nn,'_', ' ')
    nn2 <- bb_snp[Species2==nn1, Species] %>% unique()
    my_phylo$tip.label <- 
        stringr::str_replace_all(my_phylo$tip.label, nn, nn2)
}


# species_group <- unique(bb_snp$Group2) #bb_snp[,c('Species', 'Group2')] %>% unique() %>% pull(Species, Group2)
# 
# my_phylo <- groupOTU(my_phylo, species_group)

my_phylo.dt <- as_tibble(my_phylo) %>% as.data.table

plot(my_phylo)

################################################################################


# bb_snp$Group <- factor(bb_snp$Group, 
#                        levels=c('unicellular','primates', 'mammals',  'birds',  
#                                 'fish', 'worm', 'arthropods', 'fungus','plants'))

InverseTree <- inverseA(my_phylo,scale=TRUE)$Ainv
bb_snp$log_callable <- log(bb_snp$Callable)

vv <- diag(2) 
diag(vv) <- c(1e6, 1e-8)
my_prior <- list(B=list(mu=c(0,1), V=vv),
                 G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m1 <- MCMCglmm(Events ~ log_callable,
               random= ~ Species,
               data=bb_snp,
               ginverse=list(Species=InverseTree),
               prior = my_prior,
               nitt=1001000, thin=1000, burnin=10000,
               family = "poisson",
               pr=TRUE, verbose = FALSE)

HPDinterval((m1$VCV[,'Species'])/rowSums(m1$VCV))
mean((m1$VCV[,'Species'])/rowSums(m1$VCV))



# create a table for prediction
bb_snp_pred <- data.table('Species'=c(my_phylo$node.label, 
                                   my_phylo$tip.label)) %>%
    filter(Species!='n') %>%
    mutate(Events=0) %>%
    mutate(log_callable=0)

bb_snp_pred[Species=='n87', log_callable:=1]

# prediction based on model m1
bb_snp_pred <- predict(m1, interval="confidence", newdata=bb_snp_pred, 
                       type='response', marginal = NULL) %>%
    cbind(bb_snp_pred,.)

###############################################################################


plot(my_phylo)

nodelabels(text = my_phylo$node.label, cex=0.6)

## get root mean
anc0 <- bb_snp_pred[Species=='n262',4:6]
log_anc0 <- log(anc0)

bb_snp_pred <- bb_snp_pred %>% 
    mutate(t1=' [', t2=',', t3=']') %>%
    mutate(fit1 = log(fit), 
           lwr1 = log(lwr), 
           upr1 = log(upr)) %>%
    mutate(fit2 = sprintf('%.2f', fit*1e9), 
           lwr2 = sprintf('%.2f', lwr*1e9), 
           upr2 = sprintf('%.2f', upr*1e9)) %>%
    tidyr::unite('tr', fit2,t1,lwr2,t2,upr2,t3, sep = '', remove = FALSE) %>%
    tidyr::unite('CI', t1, lwr2, t2, upr2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>%
    as.data.table

tip_pred <- bb_snp_pred[!Species %like% '^n\\d+']

bb_snp_aggre <- bb_snp[,.(nn=sum(.SD$Events), NN=sum(.SD$Callable)/1e9), by=Species]

tip_pred <- left_join(tip_pred, bb_snp_aggre, by='Species') %>%
    mutate(NN2=sprintf('%.2f', NN))



scol_1 <- c('primates'=primates_col, 'mammals'=mammals_col, 'birds'=birds_col, 
            'reptiles'=reptiles_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
            'nematodes'=nematodes_col,'plants'=plants_col, 'fungi'=fungi_col,
            'others'=others_col)

#ggtree(my_phylo) + geom_text2(aes(label=node)) # show node order
node_sets <- data.table('clade' = c("primates","mammals","birds",
                                    "reptiles","fish", "arthropods",
                                    'nematodes', 'fungi',"plants"),
                        'Species' = c('n133','n144','n188',
                                      'n156','n63','n229',
                                      'n224', 'n43', 'n245'))


clade_node <-  my_phylo.dt[label %in% node_sets$Species] %>%
    rename('Species'='label') %>%
    left_join(.,node_sets, by='Species') %>%
    mutate(clade = factor(clade, levels=c("primates","mammals","birds",
                                          "reptiles","fish", "arthropods",
                                          'nematodes', 'fungi',"plants")))
   

highlight_nodes <- clade_node$node

#ggtree(my_phylo) + geom_text2(aes(label=branch.length))

f0 <- ggtree(my_phylo) + geom_tiplab(inetype='dashed', linesize=.3) +
    xlim_tree(1000) + theme_tree2() 
    
f00 <- revts(f0)


f1 <- f00 %<+%  clade_node +
    geom_hilight(mapping=aes(subset = node %in% clade_node$node, fill=clade),
                 type = "gradient") +
    scale_fill_manual(values = scol_1) +
    theme(legend.title = element_blank(), legend.position = c(0.1,0.8),
          legend.text = element_text(size=15)) 

# modify: f1 - color annotation; f00 - no color
# draw two figures and manually merge; cannot draw at once, may be a ggtree bug
    
f12 <- facet_plot(f1, panel='count', data=tip_pred,
                  geom=geom_text2,
                  mapping=aes(x=1, label=as.character(nn)),
                  hjust=1) 

f12_2 <- facet_plot(f12, panel='count', data=tip_pred,
                    geom_text2,
                    mapping=aes(x=3, label=as.character(NN2)),
                    hjust=1) + xlim_expand(c(0, 3.5), 'count') 

f123 <- facet_plot(f12_2, panel='interval', data=tip_pred,
                   geom_linerange,
                   mapping = aes(x=fit1, xmin=lwr1, xmax=upr1))

f123_2 <- facet_plot(f123, panel='interval', data=tip_pred,
                     geom_vline,
                     mapping = aes(xintercept=log_anc0$fit[1]),
                     linetype = 'dashed')

f1234 <- facet_plot(f123_2, panel='interval', data=tip_pred,
                    geom_point,
                    mapping = aes(x=fit1))

f12345 <- facet_plot(f1234, panel='mean', data=tip_pred,
                     geom_text2,
                     mapping = aes(x=1, label=fit2),
                     hjust=1) + xlim_expand(c(0,0), 'mean') 

f_all <- facet_plot(f12345, panel='ci', data=tip_pred,
                    geom_text2,
                    mapping = aes(x=1, label=CI),
                    hjust=1) + xlim_expand(c(0.85, 1.02), 'mean') 

gt = ggplot_gtable(ggplot_build(f_all))
# gtable_show_layout(gt) 
#gt # see plot layout in table format
f1_idx <- gt$layout$l[grep('panel-1-1', gt$layout$name)] 
gt$widths[f1_idx] = 1.8*gt$widths[f1_idx] 

f2_idx <- gt$layout$l[grep('panel-1-2', gt$layout$name)] 
gt$widths[f2_idx] = 0.55*gt$widths[f2_idx] 

f3_idx <- gt$layout$l[grep('panel-1-3', gt$layout$name)] 
gt$widths[f3_idx] = 0.8*gt$widths[f3_idx]

f4_idx <- gt$layout$l[grep('panel-1-4', gt$layout$name)] 
gt$widths[f4_idx] = 0.2*gt$widths[f4_idx] 

f5_idx <- gt$layout$l[grep('panel-1-5', gt$layout$name)] 
gt$widths[f5_idx] = 0.7*gt$widths[f5_idx]

pdf('tt1.pdf', width = 15, height = 20)
grid.draw(gt)
dev.off()


### get node estimates and update the figure manually

left_join(node_sets, bb_snp_pred, by='Species') %>% select(clade, fit2, lwr2, upr2)

################################################################################


####
# Add Method as a fixed term

vv <- diag(3) 
diag(vv) <- c(1e6, 1e-8,1e6)
prior_m1.2 <- list(B=list(mu=c(0,1,0), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m1.2 <- MCMCglmm(Events ~ log_callable + Method,
                 random= ~ Species,
                 data=bb_snp,
                 ginverse=list(Species=InverseTree),
                 prior = prior_m1.2,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson", verbose = FALSE,
                 pr=TRUE)

summary(m1.2)



####
# Add Reproduction as a fixed term


vv <- diag(3) 
diag(vv) <- c(1e6, 1e-8, 1e6)
prior_m1.3 <- list(B=list(mu=c(0,1,0), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m1.3 <- MCMCglmm(Events ~ log_callable + reproduction,
                 random= ~ Species,
                 data=bb_snp[reproduction %in% c('sexual','asexual')],
                 ginverse=list(Species=InverseTree),
                 prior = prior_m1.3,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson", verbose = FALSE,
                 pr=TRUE)

summary(m1.3)




# Iterations = 10001:1000001
# Thinning interval  = 1000
# Sample size  = 991
# 
# DIC: 1828.533
# 
# G-structure:  ~Species
# 
# post.mean l-95% CI u-95% CI eff.samp
# Species     3.646    2.079    5.404      991
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units    0.2285   0.1729   0.2924      991
# 
# Location effects: Events ~ log_callable + reproduction
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC
# (Intercept)         -21.8288 -23.3479 -20.3155      991 <0.001 **
#     log_callable          1.0000   0.9998   1.0002      991 <0.001 **
#     reproductionsexual    1.6397   0.9659   2.2725      991 <0.001 **
#     ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



################################################################################
################# Dmel

bb_dmel <- bb_snp %>% filter(Species == 'Drosophila melanogaster')

bb_dmel[is.na(u_lower), u_lower:=apply(.SD,1,function(x) inferLowerCI(x[1],x[2])), 
        .SDcols=c('Events','Callable')]

bb_dmel[is.na(u_upper), u_upper:=apply(.SD,1,function(x) inferUpperCI(x[1],x[2])), 
        .SDcols=c('Events','Callable')]

keep_rows = c('PID','PID2','u_mean', 'u_lower', 'u_upper', 'Events', 'Callable', 
              'Method')

bb_dmel <- bb_dmel %>% select(all_of(keep_rows)) %>% arrange(Method) %>%
    mutate(rid=10:1) %>% 
    mutate(Callable = sprintf('%.0f', Callable)) %>%
    mutate(u0=sprintf('%.2f',u_mean*1e9)) %>%
    mutate(u1=sprintf('%.2f',u_lower*1e9)) %>%
    mutate(u2=sprintf('%.2f',u_upper*1e9)) %>%
    mutate(t1='  [', t2=', ', t3=']') %>% 
    tidyr::unite(uCI, u0,  t1, u1, t2, u2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>%
    as.data.table()

## Plot Dmel 

f1 <- ggplot() + geom_text(aes(x=1,y=rid, label=PID2), data=bb_dmel, hjust=0) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=rid, label=Events), data=bb_dmel, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=rid, label=Callable), data=bb_dmel, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

ma_mean = mean(bb_dmel[Method=='MA', u_mean])
po_mean = mean(bb_dmel[Method=='PO', u_mean])

f4 <- ggplot(data=bb_dmel) + geom_point(aes(x=u_mean, y=rid)) +
    geom_vline(xintercept = c(ma_mean, po_mean),  linetype='dashed', col='red') +
    geom_linerange(aes(xmin=u_lower, xmax=u_upper, y=rid)) +
    theme_void() + 
    theme(axis.line.x = element_line(), 
          axis.text.x = element_text(), 
          axis.title.x = element_text()) +
    scale_x_continuous(name='mutation rate(1e-9)', 
                       breaks=c(2e-9,4e-9,6e-9,8e-9), 
                       labels = c(2,4,6,8))

f5 <- ggplot() + geom_text(aes(x=1,y=rid, label=uCI), data=bb_dmel, hjust=1) +
    theme_void() 

plot_grid(f1,f2,f3, f4, f5, nrow=1, axis='bt', align='hv', 
          rel_widths = c(0.4,0.2,0.3,0.3,0.4))

ggsave('dmel_MA_PO.pdf', width = 12, height = 5)

tt <- bb_dmel[,.(ee=sum(.SD[,Events]),nn=sum(as.numeric(.SD[,Callable]))),by=Method]
binconf(tt$ee,tt$nn)*1e9
# PointEst    Lower    Upper
# 5.276248 5.092053 5.467106
# 3.198268 2.374577 4.307680

######################################
################### mice

bb_mus <- bb_snp %>% filter(Species == 'Mus musculus')

bb_mus[is.na(u_lower), u_lower:=apply(.SD,1,function(x) inferLowerCI(x[1],x[2])), 
       .SDcols=c('Events','Callable')]

bb_mus[is.na(u_upper), u_upper:=apply(.SD,1,function(x) inferUpperCI(x[1],x[2])), 
       .SDcols=c('Events','Callable')]

bb_mus <- bb_mus %>% select(all_of(keep_rows)) %>% arrange(Method) %>%
    mutate(rid=5:1) %>% mutate(Callable = sprintf('%.0f', Callable)) %>%
    mutate(u0=sprintf('%.2f',u_mean*1e9)) %>%
    mutate(u1=sprintf('%.2f',u_lower*1e9)) %>%
    mutate(u2=sprintf('%.2f',u_upper*1e9)) %>%
    mutate(t1='  [', t2=', ', t3=']') %>% 
    tidyr::unite(uCI, u0,  t1, u1, t2, u2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>%
    as.data.table

## Plot mice
f1 <- ggplot() + geom_text(aes(x=1,y=rid, label=PID2), data=bb_mus, hjust=0) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=rid, label=Events), data=bb_mus, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=rid, label=Callable), data=bb_mus, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

ma_mean = mean(bb_mus[Method=='MA', u_mean])
po_mean = mean(bb_mus[Method=='PO', u_mean])

f4 <- ggplot(data=bb_mus) + geom_point(aes(x=u_mean, y=rid)) +
    geom_vline(xintercept = po_mean,  linetype='dashed', col='red') +
    geom_linerange(aes(xmin=u_lower, xmax=u_upper, y=rid)) +
    theme_void() + 
    theme(axis.line.x = element_line(), 
          axis.text.x = element_text(), 
          axis.title.x = element_text()) +
    scale_x_continuous(name='mutation rate(1e-9)', 
                       breaks=c(2e-9,4e-9,6e-9,8e-9), 
                       labels = c(2,4,6,8))

f5 <- ggplot() + geom_text(aes(x=1,y=rid, label=uCI), data=bb_mus, hjust=1) +
    theme_void() 

plot_grid(f1,f2,f3, f4, f5, nrow=1, axis='bt', align='hv', 
          rel_widths = c(0.4,0.2,0.3,0.3,0.4))

ggsave('mus_MA_PO.pdf', width = 12, height = 3)


tt <- bb_mus[,.(ee=sum(.SD[,Events]),nn=sum(as.numeric(.SD[,Callable]))),by=Method]
binconf(tt$ee,tt$nn)*1e9

# PointEst    Lower    Upper
# 5.400000 5.001696 5.830023
# 4.053192 3.807125 4.315163

################################################################################
#### Part Three
################################################################################

ff <- read_excel('mutation_rate_literature_updating3.xlsx', sheet='gs_gt') %>%
    #rename('Species2'='Species') %>%
    mutate(log_popsize = log(10^log10_popsize)) %>%
    as.data.table

cc <- left_join(bb_snp, ff, by='Species')
# cc[Species=='Homo sapiens', Group:='human']
# 
# InverseTree <- inverseA(my_phylo,scale=TRUE)$Ainv


################################################################################
### generqtion time

cc1 <- cc[Gt_Year!=0] %>% mutate(log_Gt = log(Gt_Year)) %>%
    select(Species, Events, log_callable, log_Gt, Group2, u_mean)
#cc1$Species2_phylo <- stringr::str_replace_all(cc1$Species2, ' ' ,'_')

vv <- diag(3)
diag(vv) <- c(1e6, 1e-8, 1e6)
prior_m3.1 <- list(B=list(mu=c(0,1,0), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))


m3.1 <- MCMCglmm(Events ~ log_callable + log_Gt, 
                 random= ~ Species,
                 data=cc1, 
                 ginverse=list(Species=InverseTree),
                 prior = prior_m3.1,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson", 
                 pr=TRUE, verbose = FALSE)


## Create prediction table conditional on all random levels
cc1_pred <- data.table()
for(ss in unique(cc1$Species)){
    # print(ss)
    gg <- unique(cc1[Species==ss, Group2])
    tmp <- data.table('log_Gt' = seq(-8.5,3.36,0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m3.1, newdata = tmp,  type='response', marginal = NULL) %>%
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_Gt,Species,V1) %>%
        mutate(Group2=gg)
    cc1_pred <- rbind(cc1_pred, tmp)
    
}



scol_1 <- c('primates'=primates_col, 'mammals'=mammals_col, 'birds'=birds_col, 
            'reptiles'=reptiles_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
            'nematodes'=nematodes_col,'plants'=plants_col, 'fungi'=fungi_col,
            'others'=others_col)


## log-scale y-axis
y_breaks = c(-26,-24,-22,-20,-18,-16,-14)

P_gt <- ggplot(data=cc1,aes(x=log_Gt, y=log(u_mean))) +
    geom_point(aes(col=Group2), size=3) +
    scale_x_continuous(name='Generation time (Days)',
                       breaks=c(-8,-6,-4,-2,0,2),
                       labels = c('0.1','1','7', '50','365','2700')) +
    scale_y_continuous(name='SNM rates', limits = range(y_breaks),
                       breaks = y_breaks,
                       labels=sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(values=scol_1) +
    geom_smooth(method='lm',formula = 'y ~ x', se=FALSE, linewidth=1, col='black') +
    geom_line(aes(x=log_Gt, y=log(V1), col=Group2, group=Species),
              data=cc1_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(size=13),
          axis.title = element_text(size=15))

################################################################################
###### genome size

cc2 <- cc[Gs_Mb!=0] %>% mutate(log_Gs_Mb = log(Gs_Mb))  %>%
    select(Species, Events, log_callable, log_Gs_Mb, Group2, u_mean)

vv <- diag(3)
diag(vv) <- c(1e6, 1e-8, 1e6)
prior_m3.2 <- list(B=list(mu=c(0,1,0), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m3.2 <- MCMCglmm(Events ~ log_callable + log_Gs_Mb, 
                 random= ~ Species,
                 data=cc2, 
                 ginverse=list(Species=InverseTree),
                 prior = prior_m3.2,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE, verbose = FALSE)

summary(m3.2)

cc2_pred <- data.table()

for(ss in unique(cc2$Species)){
    #print(ss)
    gg <- unique(cc2[Species==ss, Group2])
    tmp <- data.table('log_Gs_Mb' = seq(2.19, 8.04, 0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m3.2, newdata = tmp,  type='response', marginal = NULL) %>%
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_Gs_Mb,Species,V1) %>%
        mutate(Group2=gg)
    cc2_pred <- rbind(cc2_pred, tmp)
    
}


# scol_2 <- c('primates'=primates_col, 'mammals'=mammals_col,'birds'=birds_col, 
#             'reptiles'=reptiles_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
#             'nematodes'=nematodes_col, 'plants'=plants_col, 'fungus'=fungus_col,
#             'others'=others_col)

x_breaks3 = c(2,4,6,8)

P_gs <- ggplot(data=cc2, aes(x=log_Gs_Mb, y=log(u_mean))) +
    geom_point(aes(col=Group2), size=3) +
    scale_x_continuous(name='Genome size (Mb)',
                       breaks = x_breaks3,
                       labels = sprintf('%.1f', exp(x_breaks3))) +
    scale_y_continuous(name='SNM rates',
                       limits = range(y_breaks),
                       breaks = y_breaks,
                       labels = sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(name='', values=scol_1) +
    geom_smooth(method='lm', formula = 'y~x', se=FALSE, linewidth=1, col='black') +
    geom_line(aes(x=log_Gs_Mb, y=log(V1), col=Group2, group=Species),
              data=cc2_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(size=13),
          axis.title = element_text(size=15))

################################################################################
### population size


cc3 <- cc[log_popsize!=0]  %>%
    select(Species, Events, log_callable, log_popsize, Group2, u_mean)
cc3$Group2 <- factor(cc3$Group2, levels=c('primates','mammals','birds',
                                             'reptiles', 'fish','arthropods',
                                             'nematodes', 'plants', 'fungi','others'))


vv <- diag(3)
diag(vv) <- c(1e6, 1e-8, 1e6)
prior_m3.3 <- list(B=list(mu=c(0,1,0), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))


m3.3 <- MCMCglmm(Events ~ log_callable + log_popsize, 
                 random= ~ Species,
                 data=cc3, 
                 ginverse=list(Species=InverseTree),
                 prior = prior_m3.3,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE, verbose = FALSE)

summary(m3.3)

cc3_pred <- data.table()

for(ss in unique(cc3$Species)){
    #print(ss)
    gg <- unique(cc3[Species==ss, Group2])
    tmp <- data.table('log_popsize' = seq(12.66,34.95,0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m3.3, newdata = tmp,  type='response', marginal = NULL) %>%
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_popsize,Species,V1) %>%
        mutate(Group2=gg)
    cc3_pred <- rbind(cc3_pred, tmp)
}



x_breaks2 = c(15,20,25,30,35)

P_ps <- ggplot(data=cc3,aes(x=log_popsize, y=log(u_mean))) +
    geom_point(aes(col=Group2), size=3) +
    scale_x_continuous(name='Population size' ,
                       breaks = x_breaks2,
                       labels=sprintf('%.1e', exp(x_breaks2))) +
    scale_y_continuous(name='SNM rates',
                       limits = range(y_breaks),
                       breaks = y_breaks,
                       labels=sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(name='', values=scol_1, drop=FALSE) +
    geom_smooth(method='lm', formula = 'y~x', se=FALSE, linewidth=1, col='black') +
    geom_line(aes(x=log_popsize, y=log(V1), col=Group2, group=Species),
              data=cc3_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = c(0.53,0.15),
          legend.text = element_text(size=12),
          axis.text = element_text(size=13),
          axis.title = element_text(size=15)) +
    guides(color=guide_legend(nrow=3, byrow=FALSE))

################################################################################
#### all together

# all factors together

vv <- diag(5)
diag(vv) <- c(1e6, 1e-8, 1e6, 1e6, 1e6)
prior_m.ccx <- list(B=list(mu=c(0,1,0,0,0), V=vv),
                    G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

ccx <- cc[Gt_Year!=0 & log_popsize !=0] %>%
    mutate(log_Gt_Year = log(Gt_Year), log_Gs_Mb = log(Gs_Mb)) %>%
    select(Species, Events, log_callable, log_popsize, log_Gt_Year, log_Gs_Mb,
           Group2, u_mean)

m.ccx <- MCMCglmm(Events ~ log_callable + log_Gs_Mb + log_popsize + log_Gt_Year, 
                  random= ~ Species,
                  data=ccx, 
                  ginverse=list(Species=InverseTree),
                  prior = prior_m.ccx,
                  family = 'poisson',
                  nitt=1001000, thin=1000, burnin=10000,
                  pr=TRUE, verbose = FALSE)

summary(m.ccx)
## drop genome size

vv <- diag(4)
diag(vv) <- c(1e6, 1e-8, 1e6, 1e6)
prior_m.ccy <- list(B=list(mu=c(0,1,0,0), V=vv),
                    G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

ccy <- cc[Gt_Year!=0 & log_popsize !=0] %>%
    mutate(log_Gt_Year = log(Gt_Year)) %>%
    select(Species, Events, log_callable, log_popsize, log_Gt_Year,
           Group2, u_mean)

m.ccy <- MCMCglmm(Events ~ log_callable + log_popsize + log_Gt_Year, 
                  random= ~ Species,
                  data=ccx, 
                  ginverse=list(Species=InverseTree),
                  prior = prior_m.ccy,
                  family = 'poisson',
                  nitt=1001000, thin=1000, burnin=10000,
                  pr=TRUE, verbose = FALSE)

summary(m.ccy)


## ## drop census population size

vv <- diag(4)
diag(vv) <- c(1e6, 1e-8, 1e6, 1e6)
prior_m.ccz <- list(B=list(mu=c(0,1,0,0), V=vv),
                    G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))


ccz <- cc[Gt_Year!=0 & Gs_Mb!=0] %>% 
    mutate(log_Gs_Mb = log(Gs_Mb))  %>%
    mutate(log_Gt_Year = log(Gt_Year)) %>%
    select(Species, Events, log_callable, log_Gs_Mb, Group2,log_Gt_Year, u_mean)


m.ccz <- MCMCglmm(Events ~ log_callable + log_Gs_Mb + log_Gt_Year, 
                  random= ~ Species,
                  data=ccz, 
                  ginverse=list(Species=InverseTree),
                  prior = prior_m.ccz,
                  family = 'poisson',
                  nitt=1001000, thin=1000, burnin=10000,
                  pr=TRUE, verbose = FALSE)

summary(m.ccz)

################################################################################
#### nucleotide diversity
## mutation rate predication on model `m1`

mu_pred <- bb_snp_pred[!Species %like% '^n\\d+', c('fit', 'Species')]

div <- cc[diversity!=0, c('Species', 'diversity')] %>% unique()

mu_pred_div <- left_join(div, mu_pred, by = 'Species') %>%
    mutate(log_div = log(diversity)) %>%
    mutate(log_mu = log(fit)) %>%
    mutate(diversity=NULL, fit=NULL) 


# 
# ##
# # create data table for model
# d1 <- bb_snp %>% select(Species, Events, log_callable) %>%
#     mutate(trait = 'mu') %>%
#     mutate(family = 'poisson') %>%
#     rename(y = Events)
# 
# d2 <- mu_pred_div %>% select(Species, log_div) %>%
#     rename(y = log_div) %>%
#     mutate(trait = 'div') %>%
#     mutate(family = 'gaussian') %>%
#     mutate(log_callable = 0) %>%
#     select(Species, y, log_callable, trait, family)
# 
# d12 <- rbind(d1,d2)
# d12$trait <- factor(d12$trait, levels = c('mu', 'div'))
# 
# 
# ## fit bivariate model
# prior_m4.2 <- list(R=list(V=diag(2), nu=0.002),
#                    G=list(G1=list(V=diag(2), nu=2, 
#                                   alpha.V=diag(2)*100, beta.mu=0, beta.V=10)))
# 
# m4.2_biv.x <- MCMCglmm(y ~ trait - 1 + at.level(trait, 'mu'):log_callable,
#                        random = ~ante1(trait):Species,
#                        rcov = ~idh(trait):units,
#                        data = d12,
#                        ginverse = list(Species=InverseTree),
#                        prior = prior_m4.2,
#                        nitt=1001000, thin=1000, burnin=10000,
#                        family = NULL, verbose = FALSE
#                        )

# #plot(posterior.ante(m4.2_biv.x$VCV[,1:4], k=1))
# posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] %>% HPDinterval()
# 
# posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] %>% mean
# # p-value
# sum(posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] > 1) / nrow(m4.2_biv.x$VCV)
# 


###############################################################################

d1 <- bb_snp %>% select(Species, Events, Callable)
  
d1 <- d1[,.(y=sum(Events), log_callable=log(sum(Callable))),by=c('Species')] %>% 
    mutate(trait = 'mu') %>%
    mutate(family = 'poisson')

d12 <- rbind(d1,d2)
d12$trait <- factor(d12$trait, levels = c('mu', 'div'))


vv <- diag(3)
diag(vv) <- c(1e6, 1e6, 1e-8)

prior_m4.2_biv.x <- list(R=list(V=diag(2), nu=0.002, beta.mu=1, beta.V=10),
                   B=list(mu=c(0,0,1),V=vv),
                   G=list(G1=list(V=diag(2), nu=2,
                                  alpha.V=diag(2)*100, beta.mu=1, beta.V=10)))


m4.2_biv.x <- MCMCglmm(y ~ trait - 1 + at.level(trait, 'mu'):log_callable,
                       random = ~ante1(trait):Species,
                       rcov = ~ante1(trait):units,
                       data = d12,
                       ginverse = list(Species=InverseTree),
                       prior = prior_m4.2_biv.x,
                       nitt=1001000, thin=1000, burnin=10000,
                       family = NULL, verbose = FALSE
                       )

posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] %>% HPDinterval()

posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] %>% mean
# p-value
sum(posterior.ante(m4.2_biv.x$VCV[,1:4], k=1)[,3] > 1) / nrow(m4.2_biv.x$VCV)





# add a regression line
i_mu <- summary(m4.2_biv.x)$solutions[1]
i_div <- summary(m4.2_biv.x)$solutions[2]
b_slope <- posterior.ante(m4.2_biv.x$VCV[,1:4], k=1) %>% colMeans() %>% .[[3]]
b_slope2 = 1
plot_df <- data.table(plot_mu = seq(min(mu_pred_div$log_mu)-1,
                                    max(mu_pred_div$log_mu)+1,
                                    0.01),
                      plot_div = 0,
                      plot_div2 = 0)

for(i in 1:nrow(plot_df)){
    imu <- plot_df[i,plot_mu]
    y_pred <- i_div + b_slope*(imu - i_mu)
    y_pred2 <- i_div + b_slope2*(imu - i_mu)
    plot_df[i, plot_div:=y_pred]
    plot_df[i, plot_div2:=y_pred2]
}




# add group for colour aesthetic

mu_pred_div <- left_join(mu_pred_div, unique(bb_snp[,c('Species','Group2')]), by='Species')


# scol_4 <- c('primates'=primates_col, 'mammals'=mammals_col,
#             'birds'=birds_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
#             'worm'=worm_col, 'plants'=plants_col, 'fungus'=fungus_col,
#             'unicellular'=unicellular_col)


## put mutation rate on Y-axis, pi on X-axis
# mu_pred_div[Species2_phylo=='Homo_sapiens', Group:='human']

x_breaks3 = c(-25,-23,-21,-19,-17)
y_breaks3 = seq(-10,0,1)


P_div <- ggplot(data=mu_pred_div,aes(x=log_mu, y=log_div)) +
    geom_point(aes(col=Group2), size=3) +
    scale_y_continuous(name='Nucleotide diversity (%)',
                       limits = c(-8.3,-0.6),
                       breaks = y_breaks3,
                       labels=sprintf('%.2f', exp(y_breaks3)*100)) +
    scale_x_continuous(name='Predicted SNM rates',
                       limits = range(x_breaks3),
                       breaks = x_breaks3,
                       labels=sprintf('%.1e', exp(x_breaks3))) +
    geom_line(aes(x=plot_mu, y=plot_div),
              data=plot_df, linewidth=2, colour='grey', alpha=0.5) +
    geom_line(aes(x=plot_mu, y=plot_div2),
              data=plot_df, linewidth=2, colour='grey', alpha=0.5, linetype='dashed') +
    scale_color_manual(name='', values=scol_1) +
    theme_classic() +
    theme(legend.position = 'none',# c(0.85,0.4),
          axis.text = element_text(size=13),
          axis.title = element_text(size=15),
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black"))

# library(htmlwidgets)
# library(plotly)
# pp <- ggplotly(P_div, tooltip=c('x', 'y', 'Species'))
# saveWidget(pp, "p1.html")

## add generation time as a fixed effect

d12.y <- cc1[,c('Species','log_Gt')] %>% unique() %>% 
    left_join(d12, ., by='Species') %>%
    filter(!is.na(log_Gt))

vv <- diag(5)
diag(vv) <- c(1e6, 1e6, 1e6, 1e6, 1e-8)

prior_Extended_biv_model <- list(R=list(V=diag(2), nu=0.002,beta.mu=1, beta.V=10),
                         B=list(mu=c(0,0,0,0,1),V=vv),
                         G=list(G1=list(V=diag(2), nu=2,
                                        alpha.V=diag(2)*100, beta.mu=1, beta.V=10)))


Extended_biv_model <- MCMCglmm(y ~ trait - 1 + trait:log_Gt + 
                                   at.level(trait, 'mu'):log_callable ,
                               random = ~ante1(trait):Species, 
                               rcov = ~ante1(trait):units,
                               data = d12.y, 
                               ginverse = list(Species=InverseTree),
                               prior = prior_Extended_biv_model,
                               nitt=5001000, thin=1000, burnin=10000,
                               family = NULL, verbose = FALSE
                               )

#summary(Extended_biv_model)
posterior.ante(Extended_biv_model$VCV[,1:4], k=1)[,3] %>% HPDinterval()
posterior.ante(Extended_biv_model$VCV[,1:4], k=1)[,3] %>% mean
sum(posterior.ante(Extended_biv_model$VCV[,1:4], k=1)[,3] > 1) / nrow(Extended_biv_model$VCV)


#####################

no_bergeron_species <- bb_snp[!PID %like% '2023', Species] %>% unique()

d12.y.noBergeron <- d12.y[trait=='div' & Species %in% no_bergeron_species] %>% 
    rbind(d12.y[trait=='mu'],.)

Extended_biv_model.noBergeron <- MCMCglmm(y ~ trait - 1 + trait:log_Gt + 
                                   at.level(trait, 'mu'):log_callable ,
                               random = ~ante1(trait):Species, 
                               rcov = ~ante1(trait):units,
                               data = d12.y.noBergeron, 
                               ginverse = list(Species=InverseTree),
                               prior = prior_Extended_biv_model,
                               nitt=5001000, thin=1000, burnin=10000,
                               family = NULL, verbose = FALSE
                               )

posterior.ante(Extended_biv_model.noBergeron$VCV[,1:4], k=1)[,3] %>% HPDinterval()
posterior.ante(Extended_biv_model.noBergeron$VCV[,1:4], k=1)[,3] %>% mean
sum(posterior.ante(Extended_biv_model.noBergeron$VCV[,1:4], k=1)[,3] > 1) / nrow(Extended_biv_model.noBergeron$VCV)









################# some outliers?

d12.y_noUnicellular <- bb[,c('Species','Group')] %>% unique() %>% left_join(d12.y, ., by = 'Species') %>%
    filter(Group!='unicellular')


# outlier_species <- mu_pred_div[log_mu < -24, Species]
# d12.y_noOutliers <- d12.y[!(trait=='div' & Species %in% outlier_species)]
# d12.y_noOutliers <- d12.y_noOutliers[!(trait=='div' & y > -1.9)]



Extended_biv_model_noOutliers <- MCMCglmm(y ~ trait - 1 + trait:log_Gt + 
                                   at.level(trait, 'mu'):log_callable ,
                               random = ~ante1(trait):Species, 
                               rcov = ~ante1(trait):units,
                               data = d12.y_noUnicellular, 
                               ginverse = list(Species=InverseTree),
                               prior = prior_Extended_biv_model,
                               nitt=5001000, thin=1000, burnin=10000,
                               family = NULL, verbose = FALSE
                               )

posterior.ante(Extended_biv_model_noOutliers$VCV[,1:4], k=1)[,3] %>% HPDinterval()
posterior.ante(Extended_biv_model_noOutliers$VCV[,1:4], k=1)[,3] %>% mean
sum(posterior.ante(Extended_biv_model_noOutliers$VCV[,1:4], k=1)[,3] > 1) / nrow(Extended_biv_model_noOutliers$VCV)





################################################################################
### plot
plot_grid(P_gt, P_gs, P_div, P_ps, ncol=2, 
          labels = c('(A)', '(B)', '(C)', '(D)'), label_x = 0.2,
          align = 'hv')

ggsave('mu_genomeSize_generationTime_new.pdf', width = 11, height=10, units = 'in')


################################################################################
## Part Four
################################################################################

dd <- bb %>% extract(CITATION, 'Author', '^(.*?),') %>%
    tidyr::unite('PID2','Author','Date', sep=',') %>%
    filter(TYPE=='indel') %>% 
    filter(Unit=='bp/generation') %>%
    filter(!is.na(u_mean)) %>%
    mutate(u_mean = as.numeric(u_mean)) %>% 
    mutate(u_lower = as.numeric(u_lower)) %>%
    mutate(u_upper = as.numeric(u_upper)) %>%
    mutate(SE = as.numeric(SE)) %>%
    mutate(Events = as.integer(Events)) %>%
    filter(!Comments %like% 'exome') %>%
    filter(!Comments %like% 'male') %>%
    filter(!Comments %like% 'sperm') %>%
    select(c('PID', 'PID2','u_mean','u_lower','u_upper','Events','Callable', 'SE',
             'Name','Species','Method', 'Group2')) %>%
    as.data.table

# Calculate interval if SE provided
dd[is.na(Events) & is.na(u_lower) & !is.na(SE) & !is.na(u_mean), 
   u_lower:=u_mean-1.96*SE]
dd[is.na(Events) & is.na(u_upper) & !is.na(SE) & !is.na(u_mean), 
   u_upper:=u_mean+1.96*SE]

# estimate Events based on interval
dd[is.na(Events) & !is.na(u_lower), 
   Events:=apply(.SD,1,function(x) inferEvents(x[1],x[2],x[3])), 
   .SDcols=c('u_mean','u_lower','u_upper')]
# Assume 1 Event if interval no available
dd[is.na(Events) & is.na(u_lower), Events:=1]
dd[, Callable:=Events/u_mean]
dd$PID <- as.character(dd$PID)

##############################################################################


#dd[Species %in% c('Anopheles coluzzii', 'Anopheles stephensi'), Species:='Anopheles gambiae']
dd[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species:='Ostreococcus tauri']
dd[Species == 'Chlamydomonas incerta', Species:='Chlamydomonas reinhardtii']
#dd[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), Species:='Paramecium tetraurelia']
dd[Species =='Saccharomycodes ludwigii', Species:='Hanseniaspora valbyensis']


# dd[,Species2:=Species]
# dd[Species=='Amphilophus', Species2:='Amphilophus labiatus']
# dd[Species=='Anopheles stephensi', Species2:='Anopheles gambiae']
# dd[Species=='Anopheles coluzzii', Species2:='Anopheles gambiae']
# dd[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), Species2:='Paramecium tetraurelia']
# dd[Species=='Chlamydomonas incerta', Species2:='Chlamydomonas reinhardtii']
# dd[Species=='Daphnia galeata', Species2:='Daphnia dubia']
# dd[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species2:='Ostreococcus tauri']
# dd[Species=='Saccharomycodes ludwigii', Species2:='Hanseniaspora valbyensis']
# # following automatically replaced with neighbouring sppecies by treetime.org
# dd[Species=='Micromonas pusilla', Species2:='Micromonas']
# dd[Species=='Picochlorum costavermella', Species2:='Prototheca wickerhamii']
# dd[Species=='Rhodotorula toruloides', Species2:='Rhodotorula graminis']
# dd[Species=='Bombus terrestris', Species2:='Bombus']
# dd[Species=='Chironomus riparius', Species2:='Chironomus pallidivittatus']
# dd[Species=='Pristionchus pacificus', Species2:='Pristionchus']
# 
# dd_unique_species <- dd[,"Species2"] %>% unique()

################################################################################
# We extract and prune tree for indels.
# create phylogeny based on species common names




my_phylo.dd <- ape::keep.tip(my_phylo, tip=unique(dd$Species))
my_phylo.dd <- ladderize(my_phylo.dd, right = FALSE)

InverseTree.dd <- inverseA(my_phylo.dd,scale=TRUE)$Ainv

plot(my_phylo.dd)
nodelabels(text = my_phylo.dd$node.label, cex=0.6)


################################################################################
##Fit model for indels: mu ~ (1|phylogeny)

dd$log_callable <- log(dd$Callable)
vv <- diag(2) 
diag(vv) <- c(1e6, 1e-8)
prior_m5.1 <- list(B=list(mu=c(0,1), V=vv),
                   G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m5.1 <- MCMCglmm(Events ~ log_callable, 
                 random= ~ Species,
                 data=dd, 
                 ginverse=list(Species=InverseTree.dd),
                 prior = prior_m5.1,
                 nitt=1001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE, verbose = FALSE)

dd_pred <- data.table('Species'=c(my_phylo.dd$node.label, 
                               my_phylo.dd$tip.label)) %>% 
    filter(Species!='n') %>% 
    mutate(Events=0) %>%
    mutate(log_callable=ifelse(Species=='n203', 1, 0)) 

dd_pred <- predict(m5.1, marginal=NULL, interval="confidence",newdata=dd_pred) %>%
    cbind(dd_pred,.)


################################################################################

# Create data table for plot



dd_pred <- dd_pred %>% 
    mutate(t1=' [', t2=',', t3=']') %>%
    mutate(fit1 = log(fit), 
           lwr1 = log(lwr), 
           upr1 = log(upr)) %>%
    mutate(fit2 = sprintf('%.2f', fit*1e9), 
           lwr2 = sprintf('%.2f', lwr*1e9), 
           upr2 = sprintf('%.2f', upr*1e9)) %>%
    tidyr::unite('tr', fit2, t1, lwr2, t2, upr2, t3, sep = '', remove = FALSE) %>%
    tidyr::unite('CI',  t1, lwr2, t2, upr2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>%
    as.data.table

tip_pred.dd <- dd_pred[!Species %like% '^n[0-9]']

dd_aggre <- dd[,.(nn=sum(.SD$Events), NN=sum(.SD$Callable)/1e9), by=c('Species','Group2')]

tip_pred.dd <- left_join(tip_pred.dd[,-c('Events', 'log_callable')], 
                         dd_aggre, 
                         by='Species') %>%
    mutate(NN2=sprintf('%.2f', NN))

### tree mean
tree.mean.dd <- log(sum(dd$Events)/sum(dd$Callable))


# 
# dd_aggre <- dd[,.(nn=sum(.SD$Events), NN=sum(.SD$Callable)/1e9), by=Name]
# 


################################################################################

# Plot phylogentic forest tree for indels

f1 <- ggtree(my_phylo.dd) + 
    geom_tiplab(linetype='dashed', linesize=.6, offset = 0) +
    xlim_tree(2000) +
    theme_tree2() 
f1 <- revts(f1)

f12 <- facet_plot(f1, panel='count', data=tip_pred.dd,
                  geom_text2,
                  mapping=aes(x=1, label=as.character(nn)),
                  hjust=1)

f12_2 <- facet_plot(f12, panel='count', data=tip_pred.dd,
                    geom_text2,
                    mapping=aes(x=3, label=as.character(NN2)),
                    hjust=1) + xlim_expand(c(0, 3.5), 'count') 

f123 <- facet_plot(f12_2, panel='interval', data=tip_pred.dd,
                   geom_linerange,
                   mapping = aes(x=fit1, xmin=lwr1, xmax=upr1))

f123_2 <- facet_plot(f123, panel='interval', data=tip_pred.dd,
                     geom_vline,
                     mapping = aes(xintercept=tree.mean.dd),
                     linetype = 'dashed')

f1234 <- facet_plot(f123_2, panel='interval', data=tip_pred.dd,
                    geom_point,
                    mapping = aes(x=fit1))

f12345 <- facet_plot(f1234, panel='mean', data=tip_pred.dd,
                     geom_text2,
                     mapping = aes(x=1, label=fit2),
                     hjust=1) + xlim_expand(c(0,0), 'mean') 

f_all <- facet_plot(f12345, panel='ci', data=tip_pred.dd,
                    geom_text2,
                    mapping = aes(x=1, label=CI),
                    hjust=1) + xlim_expand(c(0.85, 1.02), 'mean') 


gt = ggplot_gtable(ggplot_build(f_all))
#gtable_show_layout(gt) 
#gt # see plot layout in table format
f2_idx <- gt$layout$l[grep('panel-1-2', gt$layout$name)] 
gt$widths[f2_idx] = 0.5*gt$widths[f2_idx] 

f3_idx <- gt$layout$l[grep('panel-1-3', gt$layout$name)] 
gt$widths[f3_idx] = 0.7*gt$widths[f3_idx]

f4_idx <- gt$layout$l[grep('panel-1-4', gt$layout$name)] 
gt$widths[f4_idx] = 0.2*gt$widths[f4_idx] 

f5_idx <- gt$layout$l[grep('panel-1-5', gt$layout$name)] 
gt$widths[f5_idx] = 0.5*gt$widths[f5_idx]

pdf('indel.pdf', width = 10, height = 6)
grid.draw(gt) 
dev.off()

node_sets <- data.table('clade' = c("mammals",
                                    "arthropods","fungus","plants"),
                        'Species' = c('n23','n229','n43','n245'))
left_join(node_sets, dd_pred, by='Species')

################################################################################
## Correlation Indels vs SNMs

# Use mcmcglmm to explore the correlation between INDELs and SNMs while accounting for phylogeny

dx <- bb %>% filter(TYPE %in% c('snp','indel')) %>% 
    extract(CITATION, 'Author', '^(.*?),') %>%
    tidyr::unite('PID2','Author','Date', sep=',') %>%
    filter(Unit=='bp/generation') %>%
    filter(!is.na(u_mean)) %>%
    mutate(u_mean = as.numeric(u_mean)) %>% 
    mutate(u_lower = as.numeric(u_lower)) %>%
    mutate(u_upper = as.numeric(u_upper)) %>%
    mutate(Events = as.integer(Events)) %>%
    filter(!Comments %like% 'exome') %>%
    filter(!Comments %like% 'male') %>%
    filter(!Comments %like% 'sperm') %>%
    select(c('TYPE', 'PID2','SID','u_mean','Name','Species', 'Group2')) %>%
    as.data.table


dx[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species:='Ostreococcus tauri']
dx[Species == 'Chlamydomonas incerta', Species:='Chlamydomonas reinhardtii']
dx[Species =='Saccharomycodes ludwigii', Species:='Hanseniaspora valbyensis']



dx.indel <- dx[TYPE=='indel'] %>% rename(u_indel = u_mean)
dx.snp <- dx[TYPE=='snp'] %>% rename(u_snp = u_mean)
dx_snp_indel <-  left_join(dx.indel, dx.snp[,c('SID', 'u_snp')], by='SID') %>% 
    filter(!is.na(u_snp)) %>%
    mutate(log_uindel = log(u_indel)) %>%
    mutate(log_usnp = log(u_snp)) %>%
    mutate(TYPE=NULL) 



# Prior as advised by Jarrod
prior_m5.2 <- list(
    R = list(V = diag(2), nu = 1.002), 
    G = list(
        G1= list(V = diag(2) , nu = 2, alpha.mu=c(0,0),alpha.V=diag(2)*25^2)
        
    )
)

m5.2 <- MCMCglmm(cbind(log_uindel,log_usnp) ~ trait-1, 
                 random= ~us(trait):Species,
                 rcov = ~idh(trait):units,
                 data=dx_snp_indel,
                 ginverse=list(Species=InverseTree.dd),
                 prior = prior_m5.2,
                 family = rep('gaussian',2), 
                 nitt=1001000, thin=1000, burnin=10000,
                 pr=TRUE, verbose = FALSE)

#The covariance is (significantly) positive
HPDinterval(m5.2$VCV[,2])
sum(m5.2$VCV[,2]<0)/nrow(m5.2$VCV)

HPDinterval(m5.2$VCV[,2]/(sqrt(m5.2$VCV[,1])*sqrt(m5.2$VCV[,4])))
mean(m5.2$VCV[,2]/(sqrt(m5.2$VCV[,1])*sqrt(m5.2$VCV[,4])))

#############################
# Plot correlation between indels and SNMs

XYmin = -27
XYmax = -17
xy_breaks4 = seq(XYmin,XYmax,2)

dx_snp_indel$Group2 <- factor(dx_snp_indel$Group2,
                             levels=c('primates', 'mammals', 'arthropods',
                                      'nematodes', 'fungus', 'plants', 'others'))
col4 = scol_1[levels(dx_snp_indel$Group2)]

ggplot() + geom_abline(intercept = 0, slope = 1, linetype='dashed', 
                       color='grey', size=1) +
    geom_point(data=dx_snp_indel, aes(x=log(u_snp), y=log(u_indel), color=Group2), 
               size=5) +
    scale_x_continuous(name='SNM rates', limits = c(XYmin, XYmax), 
                       breaks = xy_breaks4, label=sprintf('%.1e',exp(xy_breaks4))) +
    scale_y_continuous(name='indel mutation rates', limits = c(XYmin, XYmax), 
                       breaks = xy_breaks4, label=sprintf('%.1e',exp(xy_breaks4))) + 
    scale_color_manual(name='', values=col4) + 
    theme_classic() +
    theme(legend.position = c(0.2,0.8), legend.text = element_text(size=18),
          axis.title = element_text(size=18),
          axis.text = element_text(size=17)
    )

ggsave('snp_indel_dots.pdf', width = 7, height = 6, units = 'in')

################################################################################
###### ## Deletion and insertion


di <- bb %>% extract(CITATION, 'Author', '^(.*?),') %>%
    tidyr::unite('PID2','Author','Date', sep=',') %>%
    filter(TYPE %in% c('insertion', 'deletion')) %>% 
    filter(Unit=='bp/generation') %>% 
    filter(!is.na(u_mean)) %>%
    mutate(u_mean = as.numeric(u_mean)) %>% 
    mutate(u_lower = as.numeric(u_lower)) %>%
    mutate(u_upper = as.numeric(u_upper)) %>%
    mutate(Events = as.integer(Events)) %>%
    filter(!Comments %like% 'exome') %>%
    filter(!Comments %like% 'male') %>%
    filter(!Comments %like% 'sperm') %>%
    select(c('PID', 'PID2','TYPE','SID','u_mean','u_lower','u_upper','Events',
             'Callable','Species', 'Group2')) %>%
    as.data.table


di[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species:='Ostreococcus tauri']
di[Species == 'Chlamydomonas incerta', Species:='Chlamydomonas reinhardtii']
di[Species =='Saccharomycodes ludwigii', Species:='Hanseniaspora valbyensis']




## Assume 1 indel event where is 0
di <- di %>% mutate(Events=ifelse(Events==0,1,Events)) %>%
    mutate(Callable = Events/u_mean)


di.1 <- di %>% select('PID2', 'SID', 'Species','Group2', 'TYPE', 'Events') %>% 
    tidyr::spread(TYPE, Events) 
di.2 <- di %>% select('PID2', 'SID', 'Species','Group2', 'TYPE', 'Callable') %>% 
    tidyr::spread(TYPE, Callable) %>%
    mutate(deletion_callable = ifelse(is.infinite(deletion), insertion, deletion)) %>%
    mutate(insertion_callable = ifelse(is.infinite(insertion), deletion, insertion)) %>%
    select('SID','deletion_callable','insertion_callable')

di.12 <- left_join(di.1, di.2, by='SID') %>% as.data.table()

di.12 <- apply(di.12[,c('deletion','insertion','deletion_callable','insertion_callable')],
               1,
               function(x) DescTools::BinomRatioCI(x[1],x[3],x[2],x[4])) %>% 
    t() %>% cbind(di.12, .) %>%
    rename(c('DIrr'='V1', 'DIlr'='V2', 'DIur'='V3'))
# log transformed ratio
di.12[,DIrr:=log(DIrr)]
di.12[,DIlr:=log(DIlr)]
di.12[,DIur:=log(DIur)]

di.12[,DIrr2:=sprintf('%.2f',DIrr)]
di.12[,DIlr2:=sprintf('%.2f',DIlr)]
di.12[,DIur2:=sprintf('%.2f',DIur)]


di.12$Group2 <- factor(di.12$Group2,
                              levels=c('primates', 'mammals', 'arthropods',
                                       'nematodes', 'fungus', 'plants', 'others'))

di.123 <- di.12 %>% mutate(t1=' [', t2=',', t3=']') %>% 
    tidyr::unite('CI',  t1, DIlr2, t2, DIur2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>%
    arrange(Group2) %>%
    mutate(rid=nrow(.):1)


###############################
# Plot deletion and insertion ratio by `studies` 

f0 <- ggplot() + geom_text(aes(x=1,y=rid, label=PID2), data=di.123, hjust=1) +
    theme_void()

f1 <- ggplot() + geom_text(aes(x=1,y=rid, label=Species), data=di.123, hjust=1) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=rid, label=DIrr2), data=di.123, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=rid, label=CI), data=di.123, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f4 <- ggplot(data=di.123) + geom_point(aes(x=DIrr, y=rid)) +
    geom_vline(xintercept = 0,  linetype='dashed', col='red') +
    geom_linerange(aes(xmin=DIlr, xmax=DIur, y=rid)) +
    theme_void() + 
    theme(axis.line.x = element_line(), 
          axis.text.x = element_text(), 
          axis.title.x = element_text()) +
    scale_x_continuous(name='log ratio of deletion/insertion')

plot_grid(f0,f1,f2,f3, f4, nrow=1, axis='bt', align='hv', 
          rel_widths = c(0.4,0.2,0.1,0.2,0.3))

ggsave('deletion_insertion_ratio_full.pdf', width = 12, height = 12)


############################
# Plot deletion and insertion ratio by `species`.

di.12_aggre <- di.12[,.(deletion=sum(.SD[,deletion]),
                        insertion=sum(.SD[,insertion]),
                        deletion_callable=sum(.SD[,deletion_callable]),
                        insertion_callable=sum(.SD[,insertion_callable])),
                     by=.(Species,Group2),
                     .SDcols=c('deletion', 'insertion', 'deletion_callable', 
                               'insertion_callable')] %>%
    mutate(Group='log ratio and 95%CI')

di.12_aggre <- apply(di.12_aggre[,c('deletion','insertion','deletion_callable','insertion_callable')],
                     1,
                     function(x) DescTools::BinomRatioCI(x[1],x[3],x[2],x[4])) %>% 
    t() %>% as.data.table() %>%
    rename(c('DIrr'='V1', 'DIlr'='V2', 'DIur'='V3')) %>%  cbind(di.12_aggre, .) 

di.12_aggre[,DIrr1:=log(DIrr)]
di.12_aggre[,DIlr1:=log(DIlr)]
di.12_aggre[,DIur1:=log(DIur)]
di.12_aggre[,DIrr2:=sprintf('%.2f',log(DIrr))]
di.12_aggre[,DIlr2:=sprintf('%.2f',log(DIlr))]
di.12_aggre[,DIur2:=sprintf('%.2f',log(DIur))]


di.12_aggre <- di.12_aggre %>% mutate(t1=' [', t2=',', t3=']') %>% 
    tidyr::unite('CI',  t1, DIlr2, t2, DIur2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL) %>% arrange(Group2) %>%
    mutate(rid=nrow(.):1)


f1 <- ggplot() + geom_text(aes(x=1,y=rid, label=Species), data=di.12_aggre, hjust=1) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=rid, label=DIrr2), data=di.12_aggre, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=rid, label=CI), data=di.12_aggre, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f4 <- ggplot(data=di.12_aggre) + geom_point(aes(x=DIrr1, y=rid)) +
    geom_vline(xintercept = 0,  linetype='dashed', col='red') +
    geom_linerange(aes(xmin=DIlr1, xmax=DIur1, y=rid)) +
    theme_void() + 
    theme(axis.line.x = element_line(), 
          axis.text.x = element_text(), 
          axis.title.x = element_text()) +
    scale_x_continuous(name='log ratio of deletion/insertion')


plot_grid(f1,f2,f3, f4, nrow=1, axis='bt', align='hv', 
          rel_widths = c(0.3,0.2,0.3,0.3))
ggsave('deletion_insertion_ratio.pdf', width = 6, height = 4)


