# recombinatorics
library(tidyverse)
library(Biostrings)
codons <- names(GENETIC_CODE)
pairs <- expand_grid(MOM_CODON = codons,
                     DAD_CODON = codons)
nonpolar <- tibble(AMINO = c('G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W'),
                   TYPE = 'Nonpolar')
polar <- tibble(AMINO = c('S', 'T', 'Y', 'C', 'N', 'Q'),
                TYPE = 'Polar')
acidic <- tibble(AMINO = c('D', 'E'), TYPE = 'Acidic')
basic <- tibble(AMINO = c('K', 'H', 'R'), TYPE = 'Basic')
stop <- tibble(AMINO = c('*'), TYPE = 'STOP')

amino_type <- rbind(nonpolar, polar, acidic, basic, stop)



genetic_code_mom <- tibble(MOM_CODON = names(GENETIC_CODE),
                        AMINO_MOM = unname(GENETIC_CODE)) %>%
  left_join(amino_type %>% rename(AMINO = 'AMINO_MOM')) %>%
  rename(TYPE = 'TYPE_MOM')

genetic_code_dad <- tibble(DAD_CODON = names(GENETIC_CODE),
                        AMINO_DAD = unname(GENETIC_CODE)) %>%
  left_join(amino_type %>% rename(AMINO = 'AMINO_DAD')) %>%
  rename(TYPE = 'TYPE_DAD')


pairs <- pairs %>%
  inner_join(genetic_code_mom) %>%
  inner_join(genetic_code_dad)

recombine1 <- function(mom_codon, dad_codon){
  str_c(substr(mom_codon,1,2), substr(dad_codon, 3,3))
}

recombine2 <- function(mom_codon, dad_codon){
  str_c(substr(mom_codon, 1,1), substr(dad_codon, 2,3))
}

genetic_code1 <- tibble(RECOMB1 = names(GENETIC_CODE),
                       AMINO1 = unname(GENETIC_CODE))
genetic_code2 <- tibble(RECOMB2 = names(GENETIC_CODE),
                        AMINO2 = unname(GENETIC_CODE))

amino_type1 <- amino_type %>%
  rename(AMINO = 'AMINO1', TYPE = 'TYPE1')

amino_type2 <- amino_type %>%
  rename(AMINO = 'AMINO2', TYPE = 'TYPE2')

recomb <- pairs %>%
  mutate(RECOMB1 = recombine1(MOM_CODON, DAD_CODON),
         RECOMB2 = recombine2(MOM_CODON, DAD_CODON)) %>%
  inner_join(genetic_code1) %>%
  inner_join(genetic_code2) %>%
  inner_join(amino_type1) %>%
  inner_join(amino_type2)

diags <- recomb %>%
  filter(AMINO_MOM == AMINO_DAD)

theme_set(theme_minimal())
ggplot(diags %>%
         group_by(AMINO_MOM,TYPE1, AMINO1) %>%
         summarize(n = n()),
                   aes(x = AMINO_MOM, y= n, color = TYPE1, fill = TYPE1))+
  geom_col(alpha = 0.1)+
  geom_text(aes(label = AMINO1), position = position_stack(vjust = 0.5))+
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')

ggplot(diags %>%
         group_by(AMINO_MOM,TYPE2, AMINO2) %>%
         summarize(n = n()),
       aes(x = AMINO_MOM, y= n, color = TYPE2, fill = TYPE2))+
  geom_col(alpha = 0.1)+
  geom_text(aes(label = AMINO2), position = position_stack(vjust = 0.5))+
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')

amino_maps1 <- recomb %>%
  group_by(TYPE_MOM, TYPE_DAD, TYPE1) %>%
  summarize(n = n())

amino_maps2 <- recomb %>%
  group_by(TYPE_MOM, TYPE_DAD, TYPE2) %>%
  summarize(n = n())

ggplot(recomb, aes(x = 0, fill = TYPE1))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  facet_grid(TYPE_MOM~TYPE_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()+
  scale_fill_brewer(palette = 'Set1')

ggplot(recomb, aes(x = 0, fill = TYPE2))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  facet_grid(TYPE_MOM~TYPE_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()+
  scale_fill_brewer(palette = 'Set1')

ggplot(amino_maps2, aes(x = TYPE_DAD, y = TYPE_MOM, fill = TYPE2))+
  geom_tile()

R <- diags %>%
  filter(AMINO_MOM == 'R')
S <- diags %>%
  filter(AMINO_MOM == 'S')

ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = factor(TYPE1)))+
  geom_tile()+
  geom_text(aes(label = AMINO1, color = AMINO1))+
  labs(title = 'Upper involution')+
  theme(axis.text.x = element_text(angle = 45))
ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = AMINO2))+
  geom_tile()+
  geom_text(aes(label = AMINO2))+
  labs(title = 'Lower involution')+
  theme(axis.text.x = element_text(angle = 45))

ggplot(recomb, aes(x = 0, fill = factor(AMINO1)))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Upper ')+
  facet_grid(AMINO_MOM~AMINO_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()

ggplot(recomb, aes(x = 0, fill = factor(AMINO2)))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Lower involution')+
  facet_grid(AMINO_MOM~AMINO_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()

ggplot(recomb, aes(x = AMINO_MOM, y = AMINO_DAD, fill = AMINO1))+
  geom_tile()+
  ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Upper involution')+
  theme(axis.text.x = element_text(angle = 45))
ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = AMINO2))+
  geom_tile()+
  geom_text(aes(label = AMINO2))+
  labs(title = 'Lower involution')+
  theme(axis.text.x = element_text(angle = 45))

ggplot(recomb, aes(x = 0, fill = factor(TYPE1)))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Upper involution')+
  facet_grid(AMINO_MOM~AMINO_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()

ggplot(recomb, aes(x = 0, fill = factor(TYPE2)))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Lower involution')+
  facet_grid(AMINO_MOM~AMINO_DAD, drop = TRUE)+
  coord_polar(theta = 'y')+
  theme_void()
