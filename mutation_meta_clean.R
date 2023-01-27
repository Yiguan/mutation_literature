# This code is to analyse the germline mutation rates in Eukaryotas in context of
# phylogeny. 
# https://www.biorxiv.org/content/10.1101/2023.01.24.525323v1
# Contact: yiguan.wang@ed.ac.uk
#          yiguan.wang@outlook.com




# setwd('C://Users/ywang120/Dropbox/iMutations/')

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

human_col       = mycol[4]
primates_col    = mycol[1]
mammals_col     = mycol[2]
birds_col       = mycol[3]
fish_col        = mycol[5]
arthropods_col  = mycol[6]
worm_col        = mycol[7]
plants_col      = mycol[8] 
fungus_col      = mycol[10]
unicellular_col = mycol[11]


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

###############################################




###############################################

aa <- read_excel('mutation_rate_literature_clean.xlsx', sheet='raw') %>% 
    as.data.table()

# selection mutation accumulation (MA) studies or parent-offspring (PO) studies
bb <- aa[Method %in% c('MA','PO')]

bb[,PID] %>% unique() %>% length() # Number of studies
bb$Species %>% unique() %>% length() # Mumber of species


###############################################
### Study summery

# Number of studies by year
b_year <- bb %>% select(c('PID', 'Method')) %>% unique() %>%
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
#####
# Number of studies by species
b_species <- bb %>% select(c('PID', 'Group')) %>% unique() %>%
    group_by(Group) %>% count()

b_species$Group <- factor(b_species$Group, 
                          levels=c('human','primates','mammals','birds',
                                   'fish','arthropods','worm','plants',
                                   'fungus','unicellular'))

b_species <- b_species %>% arrange(desc(Group)) %>%
    mutate(prop = n/sum(b_species$n) * 100) %>% as.data.frame() %>%
    mutate(ypos = cumsum(prop) - 0.5*prop) %>%
    mutate(ll = sprintf('%s\n%.1f%%', Group, prop))


scol <- c('human'=human_col,'primates'=primates_col,'mammals'=mammals_col,
          'birds'=birds_col, 'fish'=fish_col,'arthropods'=arthropods_col,
          'worm'=worm_col,'plants'=plants_col, 'fungus'=fungus_col,
          'unicellular'=unicellular_col)


# Pie chart for proportion of species
P_pie_species <- ggplot(b_species, aes(x='', y=prop, fill=Group)) +
    geom_bar(width = 1, stat='identity', color='white') +
    scale_fill_manual(values=scol) +
    coord_polar('y',start=0) +
    theme_void() +
    #geom_text(aes(y=ypos, label=Group), size=5)+ 
    theme(legend.position = 'none') +
    geom_label_repel(aes(y=ypos, label=ll), size=4, nudge_x = 0.3, 
                     min.segment.length=unit(6,'lines'),
                     color='black', box.padding = 0.5, fill='white')



p12 <- plot_grid(P_pie_species, P_bar, nrow=1, labels = c( '(A)', '(B)'), label_x = 0.1)

ggsave('study_summary.jpg', p12, width=12, height=8, units = 'in', dpi = 300)
ggsave('study_summary.pdf', p12, width=12, height=8, units = 'in', dpi = 300)

##
# b1 <- bb %>% select(c('PID', 'Group', 'Method')) %>% unique()
# b1[Group!='human', Method] %>% table
# b_year %>% as.data.table %>% filter(Year>=2010) %>% 
#     select(n) %>% sum()




###############################################################################
# SNM
###############################################################################

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
    select(c('PID','PID2','u_mean','u_lower','u_upper','Events','Callable', 'SE',
             'Name','Species','Method', 'Group')) %>%
    as.data.table

# Calculate interval if SE provided

bb_snp[is.na(Events) & is.na(u_lower) & !is.na(SE) & !is.na(u_mean), 
       u_lower:=u_mean-1.96*SE]
bb_snp[is.na(Events) & is.na(u_upper) & !is.na(SE) & !is.na(u_mean), 
       u_upper:=u_mean+1.96*SE]

# estimate Events based on interval
bb_snp[is.na(Events) & !is.na(u_lower), 
       Events:=apply(.SD,1,function(x) inferEvents(x[1],x[2],x[3])), 
       .SDcols=c('u_mean','u_lower','u_upper')]
# Assume 1 Event if interval no available (only u_mean provided)
bb_snp[is.na(Events) & is.na(u_lower), Events:=1]  
## Incorrectly mentioned Saxer, 2012 in Methods in Biorxiv paper. 
## Need to update!!!!!!!!!!!!!!!!!!!!!!
##        
#        PID                                 PID2   u_mean u_lower u_upper Events Callable SE  Name      Species Method  Group
# 1: 2010005 1000 Genomes Project Consortium,2010 1.20e-08      NA      NA     NA     <NA> NA human Homo sapiens     PO  human
# 2: 2010005 1000 Genomes Project Consortium,2010 1.00e-08      NA      NA     NA     <NA> NA human Homo sapiens     PO  human
# 3: 2017014                            Yang,2017 2.17e-08      NA      NA     NA     <NA> NA maize     Zea mays     PO plants

bb_snp[, Callable:=Events/u_mean]
bb_snp$PID <- as.character(bb_snp$PID)

#####################################################################


# repair unresolved species name by treetime.org
# https://www.onezoom.org/
# Using nearest species if not found in treetime
bb_snp[,Species2:=Species]
bb_snp[Species=='Amphilophus', Species2:='Amphilophus labiatus']
bb_snp[Species=='Anopheles stephensi', Species2:='Anopheles gambiae']
bb_snp[Species=='Anopheles coluzzii', Species2:='Anopheles gambiae']
bb_snp[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), 
       Species2:='Paramecium tetraurelia']
bb_snp[Species=='Chlamydomonas incerta', Species2:='Chlamydomonas reinhardtii']
bb_snp[Species=='Daphnia galeata', Species2:='Daphnia dubia']
bb_snp[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), 
       Species2:='Ostreococcus tauri']
bb_snp[Species=='Saccharomycodes ludwigii', Species2:='Hanseniaspora valbyensis']
# following automatically replaced with neighbouring species by treetime.org
bb_snp[Species=='Micromonas pusilla', Species2:='Micromonas']
bb_snp[Species=='Picochlorum costavermella', Species2:='Prototheca wickerhamii']
bb_snp[Species=='Rhodotorula toruloides', Species2:='Rhodotorula graminis']
bb_snp[Species=='Bombus terrestris', Species2:='Bombus']
bb_snp[Species=='Chironomus riparius', Species2:='Chironomus pallidivittatus']
bb_snp[Species=='Pristionchus pacificus', Species2:='Pristionchus']


bb_unique_species <- bb_snp[,"Species2"] %>% unique()
fwrite(bb_unique_species, 'filtered_unique_species.txt', col.names = FALSE)

## go to: https://timetree.org/
## Manually make adjust branch length to avoid length of zero by changing length 0.01
########################


########################
# create phylogeny based on species common names

my_phylo <- read.tree('filtered_unique_species.nwk')

bb_snp[Group=='human', Group:='primates']
# random choose one species in each group as a proxy
bb_snp_name <- bb_snp[,.SD[1,], by='Name', .SDcols=c('Species2')] 
bb_snp_name$Species2 <- stringr::str_replace_all(bb_snp_name$Species2, ' ', '_')

my_phylo_simple <- ape::keep.tip(my_phylo, tip=bb_snp_name$Species2)
my_phylo_simple <- ladderize(my_phylo_simple)

# rename node labels
my_phylo_simple$node.label <- 
    stringr::str_replace_all(my_phylo_simple$node.label, "'", "") %>% 
    paste0('n', .)

# rename tip labels as common names
for(nn in my_phylo_simple$tip.label){
    nn2 <- bb_snp_name[Species2==nn, Name]
    my_phylo_simple$tip.label <- 
        stringr::str_replace_all(my_phylo_simple$tip.label, nn, nn2)
}




name_order <- c('amoeba','malaria parasites', 'ciliates', 'haptophyta', 'yeast', 
                        'algae',
                'human', 'chimpanzee', 'gorilla', 'orangutan', 'green monkey', 
                        'baboons', 'rhesus macaque', 'marmoset', 'owl monkey', 
                        'mouse lemur',
                'mice','cattle','cat','wolf', 'bear','platypus',
                'flycatcher',
                'cichlid fishes', 'atlantic herring',
                'nematode',
                'water flea', 'pea aphid', 'honeybee', 'bumblebee', 'butterfly', 
                        'fruit fly', 'mosquito', 'midge',
                'mushroom',
                'arabidopsis', 'cabbage', 'peach', 'white campion', 'duckweed', 
                        'rice', 'maize')

bb_snp$Name <- factor(bb_snp$Name, levels=rev(name_order))
bb_snp <- bb_snp[order(-Name)]

bb_snp$Group <- factor(bb_snp$Group, 
                       levels=c('unicellular','primates', 'mammals',  'birds',  
                                'fish', 'worm', 'arthropods', 'fungus','plants'))


####################
## MCMCglmm models


InverseTree_simple <- inverseA(my_phylo_simple,scale=FALSE)$Ainv
bb_snp$log_callable <- log(bb_snp$Callable)

vv <- diag(2) 
diag(vv) <- c(1e6, 1e-8)
my_prior <- list(B=list(mu=c(0,1), V=vv),
                 G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m1 <- MCMCglmm(Events ~ log_callable,
               random= ~ Name,
               data=bb_snp,
               ginverse=list(Name=InverseTree_simple),
               prior = my_prior,
               nitt=5001000, thin=1000, burnin=10000,
               family = "poisson",
               pr=TRUE)


bb_snp_pred <- data.table('Name'=c(my_phylo_simple$node.label, 
                                   my_phylo_simple$tip.label)) %>%
    filter(Name!='n') %>%
    mutate(Events=0) %>%
    mutate(log_callable=0)

bb_snp_pred[Name=='n9', log_callable:=1]

bb_snp_pred <- predict(m1, interval="confidence", newdata=bb_snp_pred, 
                       type='response', marginal = NULL) %>%
    cbind(bb_snp_pred,.)

####
# Add Method as a fixed term

vv <- diag(3) 
diag(vv) <- c(1e6, 1e-8,1e6)
my_prior2 <- list(B=list(mu=c(0,1,0), V=vv),
                  G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m1.x <- MCMCglmm(Events ~ log_callable + Method,
                 random= ~ Name,
                 data=bb_snp,
                 ginverse=list(Name=InverseTree_simple),
                 prior = my_prior2,
                 nitt=5001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE)

summary(m1.x)
HPDinterval(m1.x$VCV[,'Name']/(m1.x$VCV[,'Name'] + m1.x$VCV[, 'units']))





#############################
## comparison MA and PO
# Dmel

bb_dmel <- bb_snp %>% filter(Species == 'Drosophila melanogaster')

bb_dmel[,u_lower:=apply(.SD,1,function(x) inferLowerCI(x[1],x[2])), 
        .SDcols=c('Events','Callable')]

bb_dmel[,u_upper:=apply(.SD,1,function(x) inferUpperCI(x[1],x[2])), 
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
    mutate(t1=NULL, t2=NULL, t3=NULL)


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


###
# mice

bb_mus <- bb_snp %>% filter(Species == 'Mus musculus')

bb_mus[,u_lower:=apply(.SD,1,function(x) inferLowerCI(x[1],x[2])), 
       .SDcols=c('Events','Callable')]

bb_mus[,u_upper:=apply(.SD,1,function(x) inferUpperCI(x[1],x[2])), 
       .SDcols=c('Events','Callable')]

bb_mus <- bb_mus %>% select(all_of(keep_rows)) %>% arrange(Method) %>%
    mutate(rid=4:1) %>% mutate(Callable = sprintf('%.0f', Callable)) %>%
    mutate(u0=sprintf('%.2f',u_mean*1e9)) %>%
    mutate(u1=sprintf('%.2f',u_lower*1e9)) %>%
    mutate(u2=sprintf('%.2f',u_upper*1e9)) %>%
    mutate(t1='  [', t2=', ', t3=']') %>% 
    tidyr::unite(uCI, u0,  t1, u1, t2, u2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL)

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




###############
# Plot whole phylogenetic tree

## get root mean
anc0 <- bb_snp_pred[Name=='n142',4:6]
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
    mutate(t1=NULL, t2=NULL, t3=NULL)

tip_pred <- bb_snp_pred[!Name %like% '^n[0-9]']

bb_snp_aggre <- bb_snp[,.(nn=sum(.SD$Events), NN=sum(.SD$Callable)/1e9), by=Name]

tip_pred <- left_join(tip_pred, bb_snp_aggre, by='Name') %>%
    mutate(NN2=sprintf('%.2f', NN))


f1 <- ggtree(my_phylo_simple) + 
    geom_tiplab(align=TRUE, linetype='dashed', linesize=.6, offset = 50) +
    xlim_tree(1000) +
    theme_tree2() 

f1 <- revts(f1)

f12 <- facet_plot(f1, panel='count', data=tip_pred,
                  geom_text2,
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
gtable_show_layout(gt) 
#gt # see plot layout in table format
f2_idx <- gt$layout$l[grep('panel-1-2', gt$layout$name)] 
gt$widths[f2_idx] = 0.55*gt$widths[f2_idx] 

f3_idx <- gt$layout$l[grep('panel-1-3', gt$layout$name)] 
gt$widths[f3_idx] = 0.8*gt$widths[f3_idx]

f4_idx <- gt$layout$l[grep('panel-1-4', gt$layout$name)] 
gt$widths[f4_idx] = 0.2*gt$widths[f4_idx] 

f5_idx <- gt$layout$l[grep('panel-1-5', gt$layout$name)] 
gt$widths[f5_idx] = 0.7*gt$widths[f5_idx]

pdf('tt.pdf', width = 10, height = 10)
grid.draw(gt)
dev.off()

### nodes
node_sets <- data.table('clade' = c("primates_node","mammals_node","fish_node",
                                    "arthropods_node","plants_node","unicellular_node"),
                        'Name' = c('n71','n97','n53','n89','n31','n129'))

left_join(node_sets, bb_snp_pred, by='Name')



################################################################################
## Associated Factors
################################################################################

ff <- read_excel('mutation_rate_literature_clean.xlsx', sheet='gs_gt') %>%
    rename('Species2'='Species') %>%
    mutate(log_popsize = log(10^log10_popsize)) %>%
    as.data.table

cc <- left_join(bb_snp, ff, by='Species2')
cc[Species=='Homo sapiens', Group:='human']
cc$log_callable <- log(cc$Callable)

InverseTree <- inverseA(my_phylo,scale=FALSE)$Ainv

#######################
### generation time

cc1 <- cc[Gt_Year!=0] %>% mutate(log_Gt = log(Gt_Year))
cc1$Species2 <- stringr::str_replace_all(cc1$Species2, ' ' ,'_')

vv <- diag(3)
diag(vv) <- c(1e6, 1e-8,1e6)
my_prior2 <- list(B=list(mu=c(0,1,0), V=vv),
                  G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))


m2.1 <- MCMCglmm(Events ~ log_callable + log_Gt, 
                 random= ~ Species2,
                 data=cc1, 
                 ginverse=list(Species2=InverseTree),
                 prior = my_prior2,
                 nitt=5001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE)

## Create prediction table conditional on all random levels
cc1_pred <- data.table()
for(ss in unique(cc1$Species2)){
    print(ss)
    gg <- unique(cc1[Species2==ss, Group])
    tmp <- data.table('log_Gt' = seq(-5.29,3.36,0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species2' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m2.1, newdata = tmp,  type='response', marginal = NULL) %>% 
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_Gt,Species2,V1) %>%
        mutate(Group=gg)
    cc1_pred <- rbind(cc1_pred, tmp)
    
}



scol_1 <- c('human'=human_col, 'primates'=primates_col, 'mammals'='gold',
            'birds'=birds_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
            'worm'=worm_col,'plants'=plants_col, 'fungus'=fungus_col)


## log-scale y-axis
y_breaks = c(-26,-24,-22,-20,-18,-16,-14)

P_gt <- ggplot(data=cc1,aes(x=log_Gt, y=log(u_mean))) + 
    geom_point(aes(col=Group), size=3) +
    scale_x_continuous(name='Generation time (Days)', 
                       breaks=c(-4,-2,0,2), 
                       labels = c('7','50','365','2700')) +
    scale_y_continuous(name='SNM rates', limits = range(y_breaks), 
                       breaks = y_breaks, 
                       labels=sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(values=scol_1) + 
    geom_smooth(method='lm', se=FALSE, size=1, col='black') +
    geom_line(aes(x=log_Gt, y=log(V1), col=Group, group=Species2), 
              data=cc1_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text = element_text(size=13), 
          axis.title = element_text(size=15)) 


#######################
### population size

cc2 <- cc[log_popsize!=0]

cc2$Species2 <- stringr::str_replace_all(cc2$Species2, ' ' ,'_')

m2.2 <- MCMCglmm(Events ~ log_callable + log_popsize, 
                 random= ~ Species2,
                 data=cc2, 
                 ginverse=list(Species2=InverseTree),
                 prior = my_prior2,
                 nitt=5001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE)

cc2_pred <- data.table()

for(ss in unique(cc2$Species2)){
    print(ss)
    gg <- unique(cc2[Species2==ss, Group])
    tmp <- data.table('log_popsize' = seq(12.66,34.95,0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species2' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m2.2, newdata = tmp,  type='response', marginal = NULL) %>% 
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_popsize,Species2,V1) %>%
        mutate(Group=gg)
    cc2_pred <- rbind(cc2_pred, tmp)
}



scol_2 <- c('human'=human_col, 'primates'=primates_col, 'mammals'='gold',
            'birds'=birds_col, 'arthropods'=arthropods_col, 'worm'=worm_col)

x_breaks2 = c(15,20,25,30,35)

P_ps <- ggplot(data=cc2,aes(x=log_popsize, y=log(u_mean))) + 
    geom_point(aes(col=Group), size=3) +
    scale_x_continuous(name='Population size' , 
                       breaks = x_breaks2, 
                       labels=sprintf('%.1e', exp(x_breaks2))) +
    scale_y_continuous(name='SNM rates', 
                       limits = range(y_breaks), 
                       breaks = y_breaks, 
                       labels=sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(values=scol_2) +
    geom_smooth(method='lm', se=FALSE, size=1, col='black') +
    geom_line(aes(x=log_popsize, y=log(V1), col=Group, group=Species2), 
              data=cc2_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text = element_text(size=13), 
          axis.title = element_text(size=15)) 




#######################
#### genome size

cc3 <- cc[Gs_Mb!=0] %>% mutate(log_Gs_Mb = log(Gs_Mb))

cc3$Species2 <- stringr::str_replace_all(cc3$Species2, ' ' ,'_') 

m2.3 <- MCMCglmm(Events ~ log_callable + log_Gs_Mb, 
                 random= ~ Species2,
                 data=cc3, 
                 ginverse=list(Species2=InverseTree),
                 prior = my_prior2,
                 nitt=5001000, thin=1000, burnin=10000,
                 family = "poisson",
                 pr=TRUE)

cc3_pred <- data.table()

for(ss in unique(cc3$Species2)){
    print(ss)
    gg <- unique(cc3[Species2==ss, Group])
    tmp <- data.table('log_Gs_Mb' = seq(2.19,8.04, 0.01),
                      'Events' = 0,
                      'log_callable' = 0,
                      'Species2' = ss)
    tmp[1,log_callable:=1]
    tmp <- predict(m2.3, newdata = tmp,  type='response', marginal = NULL) %>% 
        cbind(tmp,.) %>% filter(log_callable==0) %>%
        select(log_Gs_Mb,Species2,V1) %>%
        mutate(Group=gg)
    cc3_pred <- rbind(cc3_pred, tmp)
    
}


scol_3 <- c('human'=human_col, 'primates'=primates_col, 'mammals'='gold',
            'birds'=birds_col, 'fish'=fish_col, 'arthropods'=arthropods_col,
            'worm'=worm_col, 'plants'=plants_col, 'fungus'=fungus_col,
            'unicellular'=unicellular_col)

x_breaks3 = c(2,4,6,8)

P_gs <- ggplot(data=cc3, aes(x=log(Gs_Mb), y=log(u_mean))) + 
    geom_point(aes(col=Group), size=3) +
    scale_x_continuous(name='Genome size (Mb)', 
                       breaks = x_breaks3, 
                       labels = sprintf('%.1f', exp(x_breaks3))) +
    scale_y_continuous(name='SNM rates', 
                       limits = range(y_breaks), 
                       breaks = y_breaks, 
                       labels = sprintf('%.1e', exp(y_breaks))) +
    scale_color_manual(name='', values=scol_3) + 
    geom_smooth(method='lm', se=FALSE, size=1, col='black') +
    geom_line(aes(x=log_Gs_Mb, y=log(V1), col=Group, group=Species2), 
              data=cc3_pred, alpha=0.3) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text = element_text(size=13), 
          axis.title = element_text(size=15))



## Plot legend 

P_legend <- ggplot(data=cc3,aes(x=log_Gs_Mb, y=log(u_mean))) + 
    geom_point(aes(col=Group), size=4) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) + 
    scale_color_manual(name='', values=scol_3) +
    theme_void() +
    theme(legend.position = c(0.2,0.5), 
          legend.text = element_text(size=15))

plot_grid(P_gt,P_ps, P_gs, P_legend, ncol=2, 
          labels = c('(A)', '(B)', '(C)'), label_x = 0.2)

ggsave('mu_genomeSize_generationTime.pdf', width = 10, height=9, units = 'in')



#######################
# all factors together

vv <- diag(5)
diag(vv) <- c(1e6, 1e-8,1e6,1e6,1e6)
my_prior3 <- list(B=list(mu=c(0,1,0,0,0), V=vv),
                  G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

cc4 <- cc[Gt_Year!=0 & log_popsize !=0] %>%
    mutate(log_Gt_Year = log(Gt_Year), log_Gs_Mb = log(Gs_Mb))

cc4$Species2 <- stringr::str_replace_all(cc4$Species2, ' ' ,'_') 

m2.4 <- MCMCglmm(Events ~ log_callable + log_Gs_Mb + log_popsize + log_Gt_Year, 
                 random= ~ Species2,
                 data=cc4, 
                 ginverse=list(Species2=InverseTree),
                 prior = my_prior3,
                 family = 'poisson',
                 nitt=5001000, thin=1000, burnin=10000,
                 pr=TRUE)






###############################################################################
## indels
###############################################################################


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
             'Name','Species','Method', 'Group')) %>%
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

dd$Species %>% unique() %>% length


dd[,Species2:=Species]
dd[Species=='Amphilophus', Species2:='Amphilophus labiatus']
dd[Species=='Anopheles stephensi', Species2:='Anopheles gambiae']
dd[Species=='Anopheles coluzzii', Species2:='Anopheles gambiae']
dd[Species %in% c('Paramecium sexaurelia','Paramecium biaurelia'), Species2:='Paramecium tetraurelia']
dd[Species=='Chlamydomonas incerta', Species2:='Chlamydomonas reinhardtii']
dd[Species=='Daphnia galeata', Species2:='Daphnia dubia']
dd[Species %in% c('Ostreococcus mediterraneus', 'Bathycoccus prasinos'), Species2:='Ostreococcus tauri']
dd[Species=='Saccharomycodes ludwigii', Species2:='Hanseniaspora valbyensis']
# following automatically replaced with neighbouring sppecies by treetime.org
dd[Species=='Micromonas pusilla', Species2:='Micromonas']
dd[Species=='Picochlorum costavermella', Species2:='Prototheca wickerhamii']
dd[Species=='Rhodotorula toruloides', Species2:='Rhodotorula graminis']
dd[Species=='Bombus terrestris', Species2:='Bombus']
dd[Species=='Chironomus riparius', Species2:='Chironomus pallidivittatus']
dd[Species=='Pristionchus pacificus', Species2:='Pristionchus']

dd_unique_species <- dd[,"Species2"] %>% unique()



###################
# create phylogeny based on species common names

dd[Group=='human', Group:='primates']
dd_name <- dd[,.SD[1,], by='Name', .SDcols=c('Species2')] 
dd_name$Species2 <- stringr::str_replace_all(dd_name$Species2, ' ', '_')

my_phylo_simple.dd <- ape::keep.tip(my_phylo, tip=dd_name$Species2)

my_phylo_simple.dd$node.label <- 
    stringr::str_replace_all(my_phylo_simple.dd$node.label, "'", "") %>% 
    paste0('n', .)

my_phylo_simple.dd <- ladderize(my_phylo_simple.dd, right = FALSE )


for(nn in my_phylo_simple.dd$tip.label){
    nn2 <- dd_name[Species2==nn, Name]
    my_phylo_simple.dd$tip.label <- 
        stringr::str_replace_all(my_phylo_simple.dd$tip.label, nn, nn2)
}


InverseTree_simple.dd <- inverseA(my_phylo_simple.dd,scale=FALSE)$Ainv
dd$log_callable <- log(dd$Callable)

m3 <- MCMCglmm(Events ~ log_callable, 
               random= ~ Name,
               data=dd, 
               ginverse=list(Name=InverseTree_simple.dd),
               prior = my_prior,
               nitt=5001000, thin=1000, burnin=10000,
               family = "poisson",
               pr=TRUE)

dd_pred <- data.table('Name'=c(my_phylo_simple.dd$node.label, 
                               my_phylo_simple.dd$tip.label)) %>% 
    filter(Name!='n') %>% 
    mutate(Events=0) %>%
    mutate(log_callable=ifelse(Name=='n87', 1, 0)) 

dd_pred <- predict(m3, marginal=NULL, interval="confidence",newdata=dd_pred) %>%
    cbind(dd_pred,.)


### get ancestor node 
anc0.dd <- dd_pred[Name=='n131',4:6]
log_anc0.dd <- log(anc0.dd)

dd_pred <- dd_pred %>% 
    mutate(t1=' [', t2=',', t3=']') %>%
    mutate(fit1 = log(fit), 
           lwr1 = log(lwr), 
           upr1 = log(upr)) %>%
    mutate(fit2 = sprintf('%.2f', fit*1e9), 
           lwr2 = sprintf('%.2f', lwr*1e9), 
           upr2 = sprintf('%.2f', upr*1e9)) %>%
    tidyr::unite('tr', fit2, t1, lwr2, t2, upr2,t3, sep = '', remove = FALSE) %>%
    tidyr::unite('CI',  t1, lwr2, t2, upr2,t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL)

tip_pred.dd <- dd_pred[!Name %like% '^n[0-9]']

dd_aggre <- dd[,.(nn=sum(.SD$Events), NN=sum(.SD$Callable)/1e9), by=Name]

tip_pred.dd <- left_join(tip_pred.dd, dd_aggre, by='Name') %>%
    mutate(NN2=sprintf('%.2f', NN))

## Plot
f1 <- ggtree(my_phylo_simple.dd) + 
    geom_tiplab(align=TRUE, linetype='dashed', linesize=.6, offset = 0) +
    xlim_tree(1000) +
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
                     mapping = aes(xintercept=log_anc0.dd$fit[1]),
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
gtable_show_layout(gt) 
#gt # see plot layout in table format
f2_idx <- gt$layout$l[grep('panel-1-2', gt$layout$name)] 
gt$widths[f2_idx] = 0.55*gt$widths[f2_idx] 

f3_idx <- gt$layout$l[grep('panel-1-3', gt$layout$name)] 
gt$widths[f3_idx] = 0.8*gt$widths[f3_idx]

f4_idx <- gt$layout$l[grep('panel-1-4', gt$layout$name)] 
gt$widths[f4_idx] = 0.2*gt$widths[f4_idx] 

f5_idx <- gt$layout$l[grep('panel-1-5', gt$layout$name)] 
gt$widths[f5_idx] = 0.7*gt$widths[f5_idx]

pdf('tt.pdf', width = 9, height = 6)
grid.draw(gt) 
dev.off()

node_sets <- data.table('clade' = c("mammals_node",
                                    "arthropods_node","plants_node"),
                        'Name' = c('n101','n89','n31'))
left_join(node_sets, dd_pred, by='Name')


###################################################
# Correlation indel and SNM

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
    select(c('TYPE', 'PID2','SID','u_mean','Name','Species', 'Group')) %>%
    as.data.table

dx.indel <- dx[TYPE=='indel'] %>% rename(u_indel = u_mean)
dx.snp <- dx[TYPE=='snp'] %>% rename(u_snp = u_mean)
dx_snp_indel <-  left_join(dx.indel, dx.snp[,c('SID', 'u_snp')], by='SID') %>% 
    left_join(., unique(bb_snp[,c('Species','Species2')]), by='Species') %>%
    filter(!is.na(u_snp)) %>%
    mutate(log_uindel = log(u_indel)) %>%
    mutate(log_usnp = log(u_snp)) %>%
    mutate(TYPE=NULL)


# get phylogenetic tree
dx_snp_indel[Group=='human', Group:='primates']
dd_name <- dx_snp_indel[,.SD[1,], by='Name', .SDcols=c('Species2')] 
dd_name$Species2 <- stringr::str_replace_all(dd_name$Species2, ' ', '_')

my_phylo_simple.dx <- ape::keep.tip(my_phylo, tip=dd_name$Species2)

my_phylo_simple.dx$node.label <- 
    stringr::str_replace_all(my_phylo_simple.dx$node.label, "'", "") %>% 
    paste0('n', .)

my_phylo_simple.dx <- ladderize(my_phylo_simple.dx, right = FALSE )

for(nn in my_phylo_simple.dx$tip.label){
    nn2 <- dd_name[Species2==nn, Name]
    my_phylo_simple.dx$tip.label <- 
        stringr::str_replace_all(my_phylo_simple.dx$tip.label, nn, nn2)
}

InverseTree_simple.dx <- inverseA(my_phylo_simple.dx,scale=FALSE)$Ainv

# Prior as advised by Jarrod
prior_d <- list(
    R = list(V = diag(2), nu = 1.002), 
    G = list(
        G1= list(V = diag(2) , nu = 2, alpha.mu=c(0,0),alpha.V=diag(2)*25^2)
        
    )
)

m4 <- MCMCglmm(cbind(log_uindel,log_usnp) ~ trait-1, 
               random= ~us(trait):Name,
               rcov = ~idh(trait):units,
               data=dx_snp_indel,
               ginverse=list(Name=InverseTree_simple.dx),
               prior = prior_d,
               family = rep('gaussian',2), 
               nitt=5001000, thin=1000, burnin=10000,
               pr=TRUE)

#The covariance is (significantly) positive
HPDinterval(m4$VCV[,2])
sum(m4$VCV[,2]<0)/nrow(m4$VCV)

HPDinterval(m4$VCV[,2]/(sqrt(m4$VCV[,1])*sqrt(m4$VCV[,4])))
mean(m4$VCV[,2]/(sqrt(m4$VCV[,1])*sqrt(m4$VCV[,4])))




XYmin = -26
XYmax = -18
xy_breaks4 = seq(XYmin,XYmax,2)

col4 = scol_3[unique(dx_snp_indel$Group)]

ggplot() + geom_abline(intercept = 0, slope = 1, linetype='dashed', 
                       color='grey', size=1) +
    geom_point(data=dx_snp_indel, aes(x=log(u_snp), y=log(u_indel), color=Group), 
               size=5) +
    scale_x_continuous(name='SNM rate', limits = c(XYmin, XYmax), 
                       breaks = xy_breaks4, label=sprintf('%.1e',exp(xy_breaks4))) +
    scale_y_continuous(name='indel mutation rate', limits = c(XYmin, XYmax), 
                       breaks = xy_breaks4, label=sprintf('%.1e',exp(xy_breaks4))) + 
    scale_color_manual(name='', values=col4) + 
    theme_classic() +
    theme(legend.position = c(0.2,0.8), legend.text = element_text(size=18),
          axis.title = element_text(size=18),
          axis.text = element_text(size=17)
    )

ggsave('snp_indel_dots.pdf', width = 7, height = 6, units = 'in')

# cor.test(dx_snp_indel$u_indel, dx_snp_indel$u_snp)




########################################################################
## insertion and deletion

d1 <- bb %>% extract(CITATION, 'Author', '^(.*?),') %>%
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
             'Callable','Name','Species','Method', 'Group')) %>%
    as.data.table


## Assume 1 indel event where is 0
d1 <- d1 %>% mutate(Events=ifelse(Events==0,1,Events)) %>%
    mutate(Callable = Events/u_mean)


d1.1 <- d1 %>% select('PID2', 'SID', 'Name','Species','Group', 'TYPE', 'Events') %>% 
    tidyr::spread(TYPE, Events) 
d1.2 <- d1 %>% select('PID2', 'SID', 'Name','Species','Group', 'TYPE', 'Callable') %>% 
    tidyr::spread(TYPE, Callable) %>%
    mutate(deletion_callable = ifelse(is.infinite(deletion), insertion, deletion)) %>%
    mutate(insertion_callable = ifelse(is.infinite(insertion), deletion, insertion)) %>%
    select('SID','deletion_callable','insertion_callable')

d1.12 <- left_join(d1.1, d1.2, by='SID') 

d12 <- apply(d1.12[,c('deletion','insertion','deletion_callable','insertion_callable')],
             1,
             function(x) DescTools::BinomRatioCI(x[1],x[3],x[2],x[4])) %>% 
    t() %>% cbind(d1.12, .) %>%
    rename(c('DIrr'='V1', 'DIlr'='V2', 'DIur'='V3'))
# log transformed ratio
d12[,DIrr:=log(DIrr)]
d12[,DIlr:=log(DIlr)]
d12[,DIur:=log(DIur)]

d12[,DIrr2:=sprintf('%.2f',DIrr)]
d12[,DIlr2:=sprintf('%.2f',DIlr)]
d12[,DIur2:=sprintf('%.2f',DIur)]



zz <- c("midge","fruit fly","bumblebee","honeybee",
        "water flea","nematode",
        "human","mice", "cattle",
        "yeast", "peach","arabidopsis","rice",
        "algae", "amoeba")

d123 <- d12 %>% mutate(t1=' [', t2=',', t3=']') %>% 
    tidyr::unite('CI',  t1, DIlr2, t2, DIur2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL, Group=NULL) %>%
    mutate(Name=factor(Name, levels=zz)) %>%
    arrange(Name) %>%
    mutate(rid=1:nrow(.))



f0 <- ggplot() + geom_text(aes(x=1,y=rid, label=PID2), data=d123, hjust=1) +
    theme_void()

f1 <- ggplot() + geom_text(aes(x=1,y=rid, label=Name), data=d123, hjust=1) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=rid, label=DIrr2), data=d123, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=rid, label=CI), data=d123, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f4 <- ggplot(data=d123) + geom_point(aes(x=DIrr, y=rid)) +
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

### aggregated deletion/insertion rates


d12_aggre <- d12[,.(deletion=sum(.SD[,deletion]),
                    insertion=sum(.SD[,insertion]),
                    deletion_callable=sum(.SD[,deletion_callable]),
                    insertion_callable=sum(.SD[,insertion_callable])),
                 by=Name,
                 .SDcols=c('Name', 'deletion', 'insertion', 'deletion_callable', 
                           'insertion_callable')] %>%
    mutate(Group='log ratio and 95%CI')

d12_aggre <- apply(d12_aggre[,c('deletion','insertion','deletion_callable','insertion_callable')],
                   1,
                   function(x) DescTools::BinomRatioCI(x[1],x[3],x[2],x[4])) %>% 
    t() %>% as.data.table() %>%
    rename(c('DIrr'='V1', 'DIlr'='V2', 'DIur'='V3')) %>%  cbind(d12_aggre, .) 

d12_aggre[,DIrr1:=log(DIrr)]
d12_aggre[,DIlr1:=log(DIlr)]
d12_aggre[,DIur1:=log(DIur)]
d12_aggre[,DIrr2:=sprintf('%.2f',log(DIrr))]
d12_aggre[,DIlr2:=sprintf('%.2f',log(DIlr))]
d12_aggre[,DIur2:=sprintf('%.2f',log(DIur))]


d12_aggre <- d12_aggre %>% mutate(t1=' [', t2=',', t3=']') %>% 
    tidyr::unite('CI',  t1, DIlr2, t2, DIur2, t3, sep = '', remove = FALSE) %>%
    mutate(t1=NULL, t2=NULL, t3=NULL, Group=NULL)

d12_aggre$Name <- factor(d12_aggre$Name, 
                         levels = rev(c("midge","fruit fly","bumblebee","honeybee",
                                         "water flea","nematode",
                                          "human","mice", "cattle",
                                          "yeast", "peach","arabidopsis","rice",
                                          "algae", "amoeba")))


f1 <- ggplot() + geom_text(aes(x=1,y=Name, label=Name), data=d12_aggre, hjust=1) +
    theme_void()

f2 <- ggplot() + geom_text(aes(x=1,y=Name, label=DIrr2), data=d12_aggre, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f3 <- ggplot() + geom_text(aes(x=1,y=Name, label=CI), data=d12_aggre, hjust=1) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))

f4 <- ggplot(data=d12_aggre) + geom_point(aes(x=DIrr1, y=Name)) +
    geom_vline(xintercept = 0,  linetype='dashed', col='red') +
    geom_linerange(aes(xmin=DIlr1, xmax=DIur1, y=Name)) +
    theme_void() + 
    theme(axis.line.x = element_line(), 
          axis.text.x = element_text(), 
          axis.title.x = element_text()) +
    scale_x_continuous(name='log ratio of deletion/insertion')


plot_grid(f1,f2,f3, f4, nrow=1, axis='bt', align='hv', 
          rel_widths = c(0.3,0.2,0.3,0.3))
ggsave('deletion_insertion_ratio.pdf', width = 6, height = 4)


