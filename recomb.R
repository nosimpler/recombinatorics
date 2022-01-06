# recombinatorics
library(tidyverse)
library(Biostrings)
alphabet <- c('A', 'T', 'G', 'C')
codons <- names(GENETIC_CODE)
pairs <- expand_grid(MOM_CODON = codons,
                     DAD_CODON = codons)

genetic_code_mom <- tibble(MOM_CODON = names(GENETIC_CODE),
                        AMINO_MOM = unname(GENETIC_CODE))
genetic_code_dad <- tibble(DAD_CODON = names(GENETIC_CODE),
                        AMINO_DAD = unname(GENETIC_CODE))

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

recomb <- pairs %>%
  mutate(RECOMB1 = recombine1(MOM_CODON, DAD_CODON),
         RECOMB2 = recombine2(MOM_CODON, DAD_CODON)) %>%
  inner_join(genetic_code1) %>%
  inner_join(genetic_code2)
ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = AMINO1))+
  geom_tile()+
  geom_text(aes(label = AMINO1))+
  labs(title = 'Upper product')+
  theme(axis.text.x = element_text(angle = 45))
ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = AMINO2))+
  geom_tile()+
  geom_text(aes(label = AMINO2))+
  labs(title = 'Lower Product')+
  theme(axis.text.x = element_text(angle = 45))

ggplot(recomb, aes(x = 0, fill = factor(AMINO1)))+
  geom_bar(position = 'fill', color = 'black')+
  #ggrepel::geom_text_repel(aes(label = AMINO1))+
  labs(title = 'Upper product')+
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
  labs(title = 'Upper product')+
  theme(axis.text.x = element_text(angle = 45))
ggplot(recomb, aes(x = MOM_CODON, y = DAD_CODON, fill = AMINO2))+
  geom_tile()+
  geom_text(aes(label = AMINO2))+
  labs(title = 'Lower Product')+
  theme(axis.text.x = element_text(angle = 45))
